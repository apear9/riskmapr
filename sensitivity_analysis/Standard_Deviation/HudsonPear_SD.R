### Script:  SENSITIVITY ANALYSIS
### Species: CYLIROSE (HUDSON PEAR)
### Varying: SD
### Notes:   N/A
### Author:  ALAN RYOU PEARSE

rm(list = ls())

## Working directory
## The [...] should point to the sensitivity_analysis folder

setwd("[...]/Data/HudsonPear/Riskfactors_susceptibility")

## Find raster files we want to process

files <- dir(pattern = ".tif$", recursive = TRUE)

## Rasters

temp <- c("DP1/Temperature.tif", "DP2/Temperature_CYLIROSE_DP2.tif", "DP3/Temperature_CYLIROSE_DP3.tif", "DP4/Temperature_CYLIROSE_DP4.tif")
rain <- c("DP1/RainAnnual.tif", "DP2/RainAnnual_CYLIROSE_DP2.tif", "DP3/RainAnnual_CYLIROSE_DP3.tif", "DP4/RainAnnual_CYLIROSE_DP4.tif")
vegl <- c("DP1/VegetationLandUse.tif", "DP2/VegetationLandUse_CYLIROSE_DP2.tif", "DP3/VegetationLandUse_CYLIROSE_DP3.tif", "DP4/VegetationLandUse_CYLIROSE_DP4.tif")
soil <- c("DP1/SoilTexture.tif", "DP2/SoilTexture_CYLIROSE_DP2.tif", "DP3/SoilTexture_CYLIROSE_DP3.tif", "DP4/SoilTexture_CYLIROSE_DP4.tif")
agoc <- c("DP1/Agochory_CYLIROSE_DP1.tif", "DP2/Agochory_CYLIROSE_DP2.tif", "DP3/Agochory_CYLIROSE_DP3.tif", "DP4/Agochory_CYLIROSE_DP4.tif")
hydc <- c("DP1/Hydrochory_CYLIROSE_DP1.tif", "DP2/Hydrochory_CYLIROSE_DP2.tif", "DP3/Hydrochory_CYLIROSE_DP3.tif", "DP4/Hydrochory_CYLIROSE_DP4.tif")
prgs <- c("DP1/PropagulesSupply_CYLIROSE_DP1.tif", "DP2/PropaguleSupply_CYLIROSE_DP2.tif", "DP3/PropaguleSupply_CYLIROSE_DP3.tif", "DP4/PropaguleSupply_CYLIROSE_DP4.tif")
zooc <- c("DP1/Zoochory_CYLIROSE_DP1.tif", "DP2/Zoochory_CYLIROSE_DP2.tif", "DP3/Zoochory_CYLIROSE_DP3.tif", "DP4/Zoochory_CYLIROSE_DP4.tif")

## Define detection periods

detp <- 1:4

## Source file containing validation functions
## Note: the [...] here should point to the sensitivity_analysis folder

source("[...]/sensitivity_functions.R")

## Loop through files and calculate Continuous Boyce Index
## The [...] here should point to the S4_Appendix folder
detections1<- readOGR("../Detections/DP1_201112.shp")
detections1@coords <- detections1@coords[,-3]
detections2 <- readOGR("../Detections/DP2_201213.shp")
detections2@coords <- detections2@coords[,-3]
detections3 <- readOGR("../Detections/DP3_201314.shp")
detections3@coords <- detections3@coords[,-3]
detections4 <- readOGR("../Detections/DP4_201510.shp")
detections4@coords <- detections4@coords[,-3]
detections5 <- readOGR("../Detections/DP5_201701.shp")
detections5@coords <- detections5@coords[,-3]
detections6 <- readOGR("../Detections/DP6_201708.shp")
detections6@coords <- detections6@coords[,-3]
dets_1 <- rbind(
  detections2,
  detections3
)
dets_2 <- rbind(
  detections3,
  detections4
)
dets_3 <- rbind(
  detections4,
  detections5
)
dets_4 <- rbind(
  detections5,
  detections6
)
dets_ <- list(dets_1, dets_2, dets_3, dets_4)
dets_[[1]]$V <- 1
dets_[[2]]$V <- 1
dets_[[3]]$V <- 1
dets_[[4]]$V <- 1
cbis_ <- vector("list", 4)

## Computation of goodness of fit statistics

for(i in detp){
  # Loop through detection periods
  temp_i <- temp[[i]]
  rain_i <- rain[[i]]
  vegl_i <- vegl[[i]]
  soil_i <- soil[[i]]
  agoc_i <- agoc[[i]]
  hydc_i <- hydc[[i]]
  prgs_i <- prgs[[i]]
  zooc_i <- zooc[[i]]
  # Loop through resolutions
  cbis_[[i]] <- numeric(length(1:50))
  for(j in 1:50){
    cbis_[[i]][j] <- compute_bn_cbi(
      persistence = c(temp_i, rain_i),
      establishment = c(vegl_i, soil_i),
      propagule = c(agoc_i, hydc_i, prgs_i, zooc_i), 
      detections = dets_[[i]],
      abundance_column = "V",
      persistence_wts = c(2, 2), # numeric vec
      persistence_sd = j + 5, # numeric
      establishment_wts = c(2, 1), # numeric vec
      establishment_sd = j + 5, # numeric
      propagule_wts = c(2, 3, 3, 2), # numeric vec
      propagule_sd = j + 5, # numeric
      suitability_sd = j, # numeric
      susceptibility_sd = j # numeric
    )$Spearman.cor
  }
}

## Plot cbis

library(ggplot2)
raster_data <- data.frame(
  cbis = c(
    cbis_[[1]],
    cbis_[[2]],
    cbis_[[3]],
    cbis_[[4]]
  ),
  detp = paste("Detection period:", rep(detp, each = 50)),
  resl = 1:50 # input standard deviation
)

ggplot(data = raster_data, aes(x = resl, y = cbis)) +
  geom_point(size = 2) +
  geom_smooth(se = FALSE) +
  facet_wrap(~ detp, nrow = 1) +
  labs(x = "Standard deviation at child nodes", y = "Continuous Boyce Index (CBI)") +
  xlim(0, NA) +
  theme_bw()

save.image("Hudson_Pear_SD.Rdata")
