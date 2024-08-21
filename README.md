# The following is a guide to the code and data used for data analysis, model simulation and figure generation.

## Raw and processed data analyzed via jupyter notebooks
    The code for analysis of raw timeseries bioreactor data (present in the form of CSVs for ODs, dilutions, 
    and single cell fluorescence) was developed in Python. The data can be loaded as a reactordata object. 
    All the functions used for analysis utilise this class. The directory contains two different kinds of code.
    The InBio GitRepo provides the libraries and class definitions that are necessary to load the bioreactor data.
    The data also contains csv files for colony counts and growth rates under processed data that is also analyzed
    through the jupyter notebook.
    The jupyter notebook consists of script and functions that were used to analyse the data once it was loaded.
    The code was developed in Windows 10 and verified by running in Window 7.
    Dependencies for the code to run without hiccups:
    - Python 3.5+
    - Pandas 1.2.4
    - Numpy 1.20
    - Matplotlib 3.4.2
    - Seaborn 0.11.1
    - SciPy 1.6.3
    - SciKit 0.24.2
    - tkinter 8.6
## Matlab data and simulation code
    Raw timeseries data was processed to obtain relevant information in individualized matlab files that were
    used to paramterize the stochastic differentiation model and simulate the dynamics given the light sequence
    (included as csv files).
    The matlab code is separated for integrated and plasmid strain. 
    Matlab code requires the statistics toolbox.
# Originally commited here: https://gitlab.inria.fr/InBio/Public/predicting-selection-effects
