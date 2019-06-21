### Script: AGGREGATE BN DATA
### Notes: YOU ONLY WANT TO RUN THIS ONCE BEFORE YOU RUN ANYTHING ELSE
### Author: ALAN RYOU PEARSE

## Clear workspace
rm(list = ls())

## Load raster package

library(raster)

## Working directory for Hudson Pear aggregations
## The [...] needs to point to the folder S4_Appendix, which is the 
## dataset for this paper. 
setwd("[...]/HudsonPear")

## Find files to load in
## These should all be .tif rasters
all_rasters <- dir(pattern = ".tif$", recursive = T)
num_rasters <- length(all_rasters)
## Load in files and aggregate
for(i in 1:num_rasters){
  current_nam <- all_rasters[i]
  current_ras <- raster(current_nam)
  for(j in 2:20){
    aggregt_ras <- aggregate(current_ras, j, fun = max)
    replace_nam <- gsub(".tif", paste0("_", j, ".tif"), current_nam)
    writeRaster(aggregt_ras, replace_nam, overwrite = TRUE)
  }
}

## Working directory for Mexican Bean Tree aggregations
setwd("[...]/MexicanBeanTree")

## Find files to load in
all_rasters <- dir(pattern = ".tif$", recursive = T)
num_rasters <- length(all_rasters)
## Load in files and aggregate
for(i in 1:num_rasters){
  current_nam <- all_rasters[i]
  current_ras <- raster(current_nam)
  for(j in 2:40){ #extend twice as much because we start off at half the resolution of the Hudson Pear rasters
    aggregt_ras <- aggregate(current_ras, j, fun = max)
    replace_nam <- gsub(".tif", paste0("_", j, ".tif"), current_nam)
    writeRaster(aggregt_ras, replace_nam, overwrite = TRUE)
  }
}

