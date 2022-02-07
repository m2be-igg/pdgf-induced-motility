# Computational model of chemotaxis calibrated with Bayesian Optimization

> PhysiCell model + Python Bayesian Optimization to replicate PDGF-induced chemotaxis migration experiments

`/model` - details on model implementation and how it can be compiled and run

`/optimization` - scripts for model optimization, reading simulation data and displaying simulation results

`/data_analysis` - notebooks used to clean experimental data and analyse experimental/computational results

## Running our optimization pipeline

The files provided in this repository can be used to **replicate our optimization pipeline using our model**, or they can be **adapted to create new PhysiCell-based optimization routines**. For example, other PhysiCell models could be integrated in the optimization process, or Bayesian Optimization could be replaced by other optimization techniques.

### Install the required Python libraries
The libraries needed to run our optimization pipeline are listed in `package-list.txt`. A new conda environment can be created from this list with the following command:

`conda create -n pdgf --file package-list.txt`

Alternatively, `vedo`, `seaborn` and `bayesian-optimization` can be installed through `pip`. These packages will automatically install some of the packages needed for our project (`numpy`, `scipy` and `pandas`).

### Compile the PhysiCell model

Run `make` inside the `model` folder.

### Run the framework

Run one of the optimization scripts provided in the `optimization` folder, using Python.