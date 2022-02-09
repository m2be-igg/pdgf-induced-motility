# This script creates a GIF showing the evolution of cells in a 3D simulation
from pathlib import Path

from vedo import Plotter, ProgressBar, screenshot

import animations
import sys
sys.path.insert(1, '../optimization')
import optimization

output_path = Path('distances-data/1ch-96h-control/replicate1/')
variables = ['position_x', 'position_y', 'position_z']
time_range = 4
time_step = 1

if __name__ == '__main__':
    # Initialize a plotter
    animation_plotter = Plotter(interactive=False)
    animations.configure_plotter(animation_plotter)
    # Initialize a progress bar to go through the time points
    progress_bar = ProgressBar(0, time_range, time_step)

    # Plot and save the results for each time point
    for timestep in progress_bar.range():
        # Get cell data
        cell_data = optimization.get_cell_data(timestep, output_path, variables)
        positions = [[x, y, z]
                     for x, y, z in zip(cell_data['position_x'],
                                        cell_data['position_y'],
                                        cell_data['position_z'])]
        # Plot cells
        animations.plot_cells(positions, animation_plotter)
        animation_plotter.show(interactive=False, resetcam=False)
        # Save the results in a PNG
        screenshot(f"day_{str(timestep).zfill(3)}.png")

        progress_bar.print(str(int(timestep)))

    # Save the final GIF
    animations.convert_into_gif("final_simulation.gif")
