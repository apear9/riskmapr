This readme file explains the contents of the sensitivity_analysis folder and their intended use. 

## Folders and contents

First notice there are a few folders:

1. Numerical_Values
2. Resolution
3. Standard_Deviation
4. Data

The first three folders correspond to the three sensitivity analyses which were run on our models for the Hudson Pear and Mexican Bean Tree. 

Each folder contains the scripts specifically relating to the associated sensitivity analysis.

The Data folder contains preprocessed data that can be used to reproduce the sensitivity analyses. 

These data are a subset of the data needed to reproduce the case studies (i.e. the data contain rasters for susceptibility risk factors for a smaller number of detection periods). 

Secondly, there are three loose .R files:

1. numerical_state_functions.R
2. preprocessing_data.R
3. sensitivity_functions.R

Only preprocessing_data.R is intended to be run by itself. 

It should be run FIRST, before you do anything else. It is a script that ensures all raster data have cells coded in the range [0, 100]. 

The other two files are sourced inside other scripts to provide useful functions. You can run them on their own but you will not see any effect because these scripts simply define functions.

## Instructions for running the code

To reproduce the Numerical_Values sensitivity analysis, 

1. Decide whether you want to do the sensitivity analysis for the Hudson Pear or the Mexican Bean Tree.
2. For the chosen species, open the associated [SPECIES_NAME]_Recoding.R file.
3. Adjust the setwd() commands / file paths as indicated in the script file and run this script from start to finish. It may take a few minutes.
4. For the chosen species, open the associated [SPECIES_NAME]_Numerical.R file. 
5. Adjust the setwd() commands / file paths as indicated in the script file and run the script from start to finish.

To reproduce the Resoluion sensitivity analysis, 

1. Open the aggregation.R file.
2. Adjust the setwd() commands / file paths as indicated in the script file. Run the script from start to finish. 
3. Decide whether you want to do the sensitivity analysis for the Hudson Pear or the Mexican Bean Tree.
4. For the chosen species, open the associated [SPECIES_NAME]_Resolution.R file.
5. Adjust the setwd() commands / file paths as indicated in the script file and run the script from start to finish.

To reproduce the Standard_Deviation sensitivity analysis, 

1. Decide whether you want to do the sensitivity analysis for the Hudson Pear or the Mexican Bean Tree.
2. For the chosen species, open the associated [SPECIES_NAME]_SD.R file.
3. Adjust the setwd() commands / file paths as indicated in the script file and run the script from start to finish. 

