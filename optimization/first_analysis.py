# This script runs the Bayesian Optimization pipeline with the full parameter space
from pathlib import Path
import subprocess

import numpy as np
from bayes_opt import BayesianOptimization, SequentialDomainReductionTransformer
from bayes_opt.logger import JSONLogger
from bayes_opt.event import Events

import optimization

FINAL_OUTPUT_PATH = Path('final_output')
EXPERIMENTAL_DATA_PATH = Path('../data-analysis/processed-data')
EXPERIMENTAL_STEM = '4mg_control_hist_day'
NUMBER_OF_REPLICATES = 3
DAYS = [1, 2, 3, 4]
PARAMS = {'sigma': (0.0, 6.0),
          'lateral_restriction': (0.0, 0.5), 'vertical_restriction': (0.5, 0.9),
          'forward_bias': (0.5, 1.0),
          'persistence_time': (10.0, 60.0),
          'cell_cell_adhesion_strength': (1.0, 10.0), 'cell_cell_repulsion_strength': (10.0, 100.0)}


def run_pipeline(sigma, lateral_restriction, vertical_restriction,
                 forward_bias,
                 persistence_time,
                 cell_cell_adhesion_strength,
                 cell_cell_repulsion_strength):

    # Run the simulation with new parameter values
    optimization.update_config_file(sigma, lateral_restriction, vertical_restriction,
                                    forward_bias,
                                    persistence_time,
                                    cell_cell_adhesion_strength,
                                    cell_cell_repulsion_strength)

    similarity_values = []

    for _ in range(NUMBER_OF_REPLICATES):
        subprocess.run('bash replicates.sh', shell=True)
        print("now creating DF")
        cells_df = optimization.read_results_into_df(FINAL_OUTPUT_PATH)
        print("now computing similarity")
        bc = optimization.compute_distance_histograms(cells_df,
                                                      EXPERIMENTAL_DATA_PATH,
                                                      EXPERIMENTAL_STEM,
                                                      DAYS)
        print(f'similarity: {bc}')

        optimization.remove_dir(Path('final_output'))

        similarity_values.append(bc)

    similarity = np.mean(similarity_values)

    return similarity


if __name__ == '__main__':
    #logger = JSONLogger(path="./logs.json")
    bounds_transformer = SequentialDomainReductionTransformer()
    optimizer = BayesianOptimization(f=run_pipeline,
                                     pbounds=PARAMS,
                                     verbose=2,
                                     random_state=1,
                                     bounds_transformer=bounds_transformer)

    #optimizer.subscribe(Events.OPTIMIZATION_STEP, logger)
    optimizer.maximize(init_points=3,
                       n_iter=4)
