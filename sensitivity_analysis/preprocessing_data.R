### Clean and preprocess the data
### Set invalid numbers to NA, etc.
### Run this only ONCE at the VERY BEGINNING; i.e. BEFORE RUNNING ANYTHING ELSE
### Alan R. Pearse

## Working directory
## You need to set this to the folder 
## sensitivity_analysis/Data

setwd("[...]/sensitivity_analysis/Data")

## Clear workspace

rm(list = ls())

## Find files to load in
## The regular expression ".tif$" should isolate only the .tif rasters
## and nothing else. 
## I can't guarantee this will always work, especially if you have
## done anything to the files in sensitivity_analysis/Data. 
## It worked for me, though. 

all_rasters <- dir(".", ".tif$", recursive = T, full.names = T)

# REDUNDANT BELOW
# ---------------------------------------------------
# all_rasters <- all_rasters[
#   -grep(".zip", all_rasters)
# ] # prepare to load everything EXCEPT .zip archives
# ---------------------------------------------------
# END OF REDUNDANT CODE

## Load raster package.
## Loop through files, replace any number above 100 with NA 
## Important because some rasters have 128 in cells which
## should really be NA. 

library(raster)
for(file in all_rasters){
  
  ras <- raster(file)
  ras[ras[] > 100] <- NA
  writeRaster(ras, file, overwrite = TRUE)
  
}

## Nice, it's finished. 
## You can of course run this multiple times without any ILL effect, but it will be a waste of time. 