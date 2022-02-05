# This script plots the cells in 3D inside a microfluidic chip for a given time point
from pathlib import Path

from vedo import Plotter

import animations
import sys
sys.path.insert(1, '../optimization')
import optimization

output_path = Path('dist-data/1ch-96h-control/replicate1/')
variables = ['position_x', 'position_y', 'position_z']
timestep = 4

if __name__ == '__main__':
    # Initialize a plotter
    animation_plotter = Plotter(interactive=False)
    animations.configure_plotter(animation_plotter)
    # Get cell data
    cell_data = optimization.get_cell_data(timestep, output_path, variables)
    positions = [[x, y, z]
                 for x, y, z in zip(cell_data['position_x'],
                                    cell_data['position_y'],
                                    cell_data['position_z'])]
    # Plot cells
    animations.plot_cells(positions, animation_plotter)
    animation_plotter.show(interactive=True)
