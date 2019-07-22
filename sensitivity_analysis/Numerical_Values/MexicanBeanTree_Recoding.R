### Script:  RECODING NUMERICAL STATES
### Species: CECROPIA (MEXICAN BEAN TREE)
### Varying: NUMERICAL STATES
### Notes:   N/A
### Author:  ALAN RYOU PEARSE

rm(list = ls())

## Working directory
## [...] needs to point to sensitivity_analysis (the root folder for the sensitivity analysis code) 

setwd("[...]/Data/MexicanBeanTree/Riskfactors_susceptibility")

## Source functions
## [...] needs to point to sensitivity_analysis (the root folder)

source("[...]/numerical_state_functions.R")
source("[...]/sensitivity_functions.R")

## Delete any files remaining from previous run(s)

to_delete <- dir(pattern = "Grid_Levels", recursive = T)
file.remove(to_delete)

## Find files to recode

files <- dir(pattern = ".tif", recursive = TRUE)

## Recoding loop

n_files <- length(files)
for(i in 1:n_files){
  r <- raster(files[i])
  nv <- length(unique(r))
  proposal_grid <- expand_grid(nv)
  proposal_grid <- cut_down_grid(proposal_grid)
  replace(r, proposal_grid)
}

