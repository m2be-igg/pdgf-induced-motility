# This script runs the Bayesian Optimization pipeline with some fixed values to reduce the parameter space
from pathlib import Path
import subprocess

import numpy as np
from bayes_opt import BayesianOptimization, SequentialDomainReductionTransformer
from bayes_opt.logger import JSONLogger
from bayes_opt.event import Events

import optimization

FINAL_OUTPUT_PATH = Path('final_output')
EXPERIMENTAL_DATA_PATH = Path('experimental/processed_data')
EXPERIMENTAL_STEM = '4mg_control_hist_day'
NUMBER_OF_REPLICATES = 3
DAYS = [1, 2, 3, 4]
PARAMS = {'sigma': (0.0, 6.0),
          'forward_bias': (0.5, 1.0),
          'persistence_time': (10.0, 60.0)}


def run_pipeline(sigma: float,
                 forward_bias: float,
                 persistence_time: float):
    """Runs the optimization pipeline with fixed values for some of the parameters"""

    cell_cell_adhesion_strength = 8.5
    cell_cell_repulsion_strength = 52.0
    vertical_restriction = 0.89
    lateral_restriction = 0.36

    # Run the simulation with new parameter values
    update_config_file(sigma, lateral_restriction, vertical_restriction, forward_bias,
                       persistence_time,
                       cell_cell_adhesion_strength,
                       cell_cell_repulsion_strength)

    similarity_values = []

    for _ in range(NUMBER_OF_REPLICATES)):
        subprocess.run('bash replicates.sh', shell=True)
        print("now creating DF")
        cells_df = read_results_into_df(FINAL_OUTPUT_PATH)
        print("now computing similarity")
        bc = optimization.compute_distance_histograms(cells_df,
                                                      EXPERIMENTAL_DATA_PATH,
                                                      EXPERIMENTAL_STEM,
                                                      DAYS)
        print(f'similarity: {bc}')

        remove_dir(Path('final_output'))

        similarity_values.append(bc)

    similarity = np.mean(similarity_values)

    return similarity


if __name__ == '__main__':
    logger = JSONLogger(path="./logs.json")
    bounds_transformer = SequentialDomainReductionTransformer()
    optimizer = BayesianOptimization(f=run_pipeline,
                                     pbounds=PARAMS,
                                     verbose=2,
                                     # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
                                     random_state=1,
                                     bounds_transformer=bounds_transformer)

    # Probe the domain with the previous optimal values
    optimizer.probe(
        params={'sigma': 2.8,
                'forward_bias': 0.56, 'persistence_time': 49.8},
        lazy=True
    )

    optimizer.subscribe(Events.OPTIMIZATION_STEP, logger)
    optimizer.maximize(init_points=3,
                       n_iter=25)
