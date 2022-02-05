# Experimental and computational data analysis

## Data processing notebooks

We provide the notebooks used to conduct our data analysis, both for experiments and computations.

- `computational-processing.ipynb` - A notebook to visualize computational results (as distance histograms and boxplots)

- `data-processing.ipynb` - A notebook used to clean experimental data prior to data analysis, and to convert distance data into histograms to be used for comparison with computational data;

- `intercellular-distance.ipynb` - A notebook created to compute the distance between cells, with some further statistical analysis and studies on the effect of the number of neighbors on intercellular distance

- `statistical-analysis.ipynb` - A notebook to run statistical tests to confirm the effect of the presence of PDGF on fibroblast motility

## Animation scripts

We also developed some scripts to visualize our simulated microfluidic setup, using the [vedo](https://vedo.embl.es/) package.

- `animations.py` - A module that provides functions to plot animations in 3D and save them as a .GIF file

- `visualization_gif.py` - A script to load and animate all the output files of a simulation, converting and saving them into a .GIF file;
- `visualization_single.py` - A script to load and animate a single time point of the simulation, displaying the result in a new window;
