# Python-based optimization pipeline

We created **Python optimization pipeline that enables users to calibrate their PhysiCell models with experimental data**, using Bayesian Optimization. Our approach relies on a publicly available [Bayesian Optimization package](https://github.com/fmfn/BayesianOptimization), `numpy`, `pandas` and `scipy`.

The scripts provided in this repository can be used in combination with our model, but they are not limited to this specific project.

## Optimization scripts

- `optimization.py` - A module that provides functions to:
  - **read cell data** from output files;
  - **save the cell data into a DataFrame** (currently working for the cells' y position);
  - **read experimental data** (currently works for our type of results);
  - **compute the distance histograms**, to be used to compute the similarity between experimental and computational data;
  - **update the XML config file** (currently working for our extension's parameters, `cell_cell_adhesion_strength` and `cell_cell_repulsion_strength`);
  - **clean the current directory** to remove unwanted files from previous runs;

- `first_analysis.py` - A script that runs a first analysis with no known initial optimal points and that varies a large number of parameters;
- `refined_analysis.py` - A script that runs the optimization pipeline for a smaller number of parameters and probes a region of interest at the start of the pipeline;
- `run_single.py` - A script that runs the base model of the pipeline, without going into optimization;
