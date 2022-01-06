# This module includes functions to ease the optimization process
import subprocess
from pathlib import Path
from xml.etree import ElementTree
from typing import List, Dict

import numpy as np
import pandas as pd

import physipy

NUMBER_OF_BINS = 30
HIST_RANGE = 650


def update_config_file(sigma: float,
                       lateral_restriction: float, vertical_restriction: float,
                       forward_bias: float,
                       persistence_time: float,
                       cell_cell_adhesion_strength: float, cell_cell_repulsion_strength: float) -> None:
    """Updates configuration file with the specified input values."""

    custom_values = {'sigma': sigma,
                     'lateral_restriction': lateral_restriction,
                     'vertical_restriction': vertical_restriction,
                     'forward_bias': forward_bias}

    mechanics_values = {'cell_cell_adhesion_strength': cell_cell_adhesion_strength,
                        'cell_cell_repulsion_strength': cell_cell_repulsion_strength}

    motility_values = {'persistence_time': persistence_time}

    custom_stem = 'user_parameters'
    mechanics_stem = 'cell_definitions/cell_definition[@name="default"]/phenotype/mechanics'
    motility_stem = 'cell_definitions/cell_definition[@name="default"]/phenotype/motility'

    file_path = Path('config/PhysiCell_settings.xml')
    tree = ElementTree.parse(file_path)

    for key, value in zip(custom_values.keys(), custom_values.values()):
        tree.find(f'{custom_stem}/{key}').text = str(value)

    for key, value in zip(mechanics_values.keys(), mechanics_values.values()):
        tree.find(f'{mechanics_stem}/{key}').text = str(value)

    for key, value in zip(motility_values.keys(), motility_values.values()):
        tree.find(f'{motility_stem}/{key}').text = str(value)

    tree.write(file_path)


def read_experimental_histograms(experimental_path: Path, file_stem: str,
                                 days: List[int]) -> Dict[int, np.ndarray]:
    """Loads the experimental data for the traveled distance histograms into a NumPy array."""
    experimental_histograms = {day: np.loadtxt(experimental_path / f'{file_stem}_{day}.csv')
                               for day in days}

    return experimental_histograms


def clean_environment(clean_data: bool = True, remove_compiled_files: bool = False) -> None:
    """Removes all previously compiled files and output data"""
    if clean_data:
        subprocess.run('make data-cleanup', shell=True)
    if remove_compiled_files:
        subprocess.run('make clean; make', shell=True)


def read_results_into_df(output_path: Path) -> pd.DataFrame:
    """Converts the simulation results for the y position into a DataFrame"""
    variables = ['position_y']
    all_cells = []
    replicates = list(output_path.glob('output*'))
    number_of_output_files = len(list(replicates[0].glob('output*.xml')))

    for replicate in replicates:
        for timestep in range(number_of_output_files):
            cell_data = physipy.get_cell_data(timestep, replicate, variables)
            cell_df = pd.DataFrame(cell_data)
            cell_df['timestep'] = timestep
            cell_df['replicate'] = replicate
            all_cells.append(cell_df)

    cells_df = pd.concat([df for df in all_cells])

    return cells_df


def compute_distance_histograms(cells_df: pd.DataFrame, experimental_path: Path, 
                                experimental_stem: str, days: List[int]) -> np.ndarray:
    """Computes the BC between the simulation and experimental results"""
    bins = np.linspace(0, HIST_RANGE, NUMBER_OF_BINS)
    similarity_coefficients = []
    computational_histograms = {}

    experimental_histograms = read_experimental_histograms(experimental_path, experimental_stem, days)

    for day in days:
        cells = cells_df[cells_df['timestep'] == day]
        distances = cells['position_y']
        weights = np.ones_like(distances) / len(distances)
        computational_histogram, bins = np.histogram(distances, weights=weights, bins=bins)
        computational_histograms[day] = computational_histogram
        similarity_matrix = [np.sqrt(comp_bin * exp_bin)
                             for comp_bin, exp_bin
                             in zip(computational_histogram, experimental_histograms[day])]

        bc = np.sum(similarity_matrix)
        similarity_coefficients.append(bc)

    return np.mean(similarity_coefficients)


def remove_dir(dir_path: Path) -> None:
    """Removes the passed directory"""
    for child in dir_path.glob('*'):
        if child.is_file():
            child.unlink()
        else:
            remove_dir(child)
    dir_path.rmdir()
