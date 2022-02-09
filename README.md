# Computational model of chemotaxis calibrated with Bayesian Optimization

> PhysiCell model + Python Bayesian Optimization to replicate PDGF-induced chemotaxis migration experiments

`/model` - details on model implementation and how it can be compiled and run

`/optimization` - scripts for model optimization, reading simulation data and displaying simulation results

`/data_analysis` - notebooks used to clean experimental data and analyse experimental/computational results

## Running our optimization pipeline

The files provided in this repository can be used to **replicate our optimization pipeline using our model**, or they can be **adapted to create new PhysiCell-based optimization routines**. For example, other PhysiCell models could be integrated in the optimization process, or Bayesian Optimization could be replaced by other optimization techniques.

### Install the required Python libraries
The libraries needed to run our optimization pipeline are listed in `requirements.txt`. We recommend using Python 3.7 or lower. A new `conda` environment with the following commands:

`conda create -n pdgf python=3.7 pip`

`pip install -r requirements.txt`

Alternatively, install the following libraries in your environment using `pip`: `bayesian-optimization`, `numpy`, `pandas` and `scipy`. 

To run the data analysis scripts, you will also need: `jupyter`, `matplotlib`, `seaborn`, `statsmodels` and `vedo`.

*This has been tested on Windows and Linux systems.*

### Compile the PhysiCell model

Run `make` inside the `model` folder.

### Run the framework

Run one of the optimization scripts provided in the `optimization` folder, using Python.
