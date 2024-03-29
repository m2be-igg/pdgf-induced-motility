# This script runs the full simulation pipeline for a set of parameter values
from pathlib import Path
import subprocess

import numpy as np

import optimization

FINAL_OUTPUT_PATH = Path('final_output')
EXPERIMENTAL_DATA_PATH = Path('../data-analysis/processed-data')
EXPERIMENTAL_STEM = '4mg_control_hist_day'
NUMBER_OF_REPLICATES = 3
DAYS = [1, 2, 3, 4]


def run_pipeline(sigma: float,
                 lateral_restriction: float, vertical_restriction: float,
                 forward_bias: float,
                 persistence_time: float,
                 cell_cell_adhesion_strength: float,
                 cell_cell_repulsion_strength: float):
    # Run the simulation with the given parameter values
    optimization.update_config_file(sigma, lateral_restriction, vertical_restriction, forward_bias,
                                    persistence_time,
                                    cell_cell_adhesion_strength,
                                    cell_cell_repulsion_strength)

    # Create a list to store the similarity values for each model replicate
    similarity_values = []

    for i in range(NUMBER_OF_REPLICATES):
        # Run the model
        subprocess.run('bash replicates.sh', shell=True)

        # Save the output data into a DataFrame
        print("now creating DF")
        cells_df = optimization.read_results_into_df(FINAL_OUTPUT_PATH)
        #cells_df.to_csv(f'./saved_dfs/dist_measure_{i}.csv', index=False)

        # Compute the BC between the model results and the experimental data
        print("now computing similarity")
        similarity = optimization.compute_distance_histograms(cells_df,
                                                              EXPERIMENTAL_DATA_PATH,
                                                              EXPERIMENTAL_STEM,
                                                              DAYS)
        print(f'similarity: {similarity}')
        similarity_values.append(similarity)

        optimization.remove_dir(FINAL_OUTPUT_PATH)

    # Print the results for all the replicates
    print(f'mean BC: {np.mean(similarity_values)}')
    print(f'std BC: {np.std(similarity_values)}')

    return similarity_values


if __name__ == '__main__':
    sigma = 2.8
    lateral_restriction = 0.36
    vertical_restriction = 0.89
    forward_bias = 0.56
    persistence_time = 49.8
    cell_cell__adhesion_strength = 8.5
    cell_cell__repulsion_strength = 52.0

    bc_coefficient = run_pipeline(sigma,
                                  lateral_restriction, vertical_restriction,
                                  forward_bias,
                                  persistence_time,
                                  cell_cell__adhesion_strength,
                                  cell_cell__repulsion_strength)
