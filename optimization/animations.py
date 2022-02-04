# This script plots cells inside a microfluidic chip in 3D for a given time point
import glob
import os
from typing import List
from PIL import Image

from vedo import Plotter, Spheres, Box


def configure_plotter(plotter: Plotter) -> None:
    """Creates a plotter with camera settings to show the chip from the top"""
    plotter.camera.SetPosition([0.0, 350.0, 2500.042])
    plotter.camera.SetFocalPoint([0.0, 350.0, 141.609])
    plotter.camera.SetViewUp([0.00, 0.0, 0.0])
    plotter.camera.SetDistance(2400.518)
    plotter.camera.SetClippingRange([1210.469, 2028.202])


def plot_cells(cell_positions: List[List[float]], plotter: Plotter) -> None:
    """Plots cells as spheres inside a bounded box, representing a microfluidic chip"""
    # Define the objects to be represented
    cells = Spheres(cell_positions, r=8, c='blue')
    chip = Box((0, 350, 0), width=700, length=1200, height=300,
               c='gray', alpha=0.5)

    # Add objects to the plotter and render the result
    plotter += cells
    plotter += chip


def convert_into_gif(file_name: str) -> None:
    """Converts the png files in the directory into a GIF"""
    # Create the frames
    frames = []
    snapshots = sorted(glob.glob("*.png"))
    for i in snapshots:
        new_frame = Image.open(i)
        frames.append(new_frame)
        os.remove(i)

    # Save into a GIF
    frames[0].save(file_name, format='GIF',
                   append_images=frames[1:],
                   save_all=True,
                   duration=300, loop=0)
