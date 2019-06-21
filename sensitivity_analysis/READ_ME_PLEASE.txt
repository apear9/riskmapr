This readme file explains the contents of the sensitivity_analysis folder and their intended use. 

First notice there are a few folders:

1. Numerical_Values
2. Resolution
3. Standard_Deviation

These folders correspond to the three sensitivity analyses which were run on our models for the Hudson Pear and Mexican Bean Tree. Each folder contains the scripts specifically relating to the associated sensitivity analysis.

Secondly, there are three loose .R files:

1. numerical_state_functions.R
2. preprocessing_data.R
3. sensitivity_functions.R

Only (2) is a script that is intended to be run by itself. It should be run FIRST, before you do anything else. The other two files are sourced inside other scripts to provide useful functions.

