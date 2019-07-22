### Script:  SENSITIVITY ANALYSIS
### Species: CECROPIA (MEXICAN BEAN TREE)
### Varying: SD
### Notes:   N/A
### Author:  ALAN RYOU PEARSE

rm(list = ls())

## Working directory
## The [...] should point to the sensitivity_analysis folder

setwd("[...]/Data/MexicanBeanTree/Riskfactors_susceptibility")

## Find raster files we want to process

files <- dir(pattern = ".tif$", recursive = TRUE)

## Rasters

elev <- c("DP1/Elevation.tif", "DP2/Elevation_DP2.tif", "DP3/Elevation_CECRSPPP_DP3.tif", "DP4/Elevation_CECRSPPP_DP4.tif", "DP5/Elevation_CECRSPPP_DP5.tif", "DP6/Elevation_CECRSPPP_DP6.tif")
rain <- c("DP1/RainDriestMonth.tif", "DP2/RainDriestMonth_DP2.tif", "DP3/RainDriestMonth_CECRSPPP_DP3.tif", "DP4/RainDriestMonth_CECRSPPP_DP4.tif", "DP5/RainDriestMonth_CECRSPPP_DP5.tif", "DP6/RainDriestMonth_CECRSPPP_DP6.tif")
vegl <- c("DP1/VegetationLandUse.tif", "DP2/VegetationLandUse_DP2.tif", "DP3/VegetationLandUse_CECRSPPP_DP3.tif", "DP4/VegetationLandUse_CECRSPPP_DP4.tif", "DP5/VegetationLandUse_CECRSPPP_DP5.tif", "DP6/VegetationLandUse_CECRSPPP_DP6.tif")
soil <- c("DP1/SoilMoisture.tif", "DP2/SoilMoisture_DP2.tif", "DP3/SoilMoisture_CECRSPPP_DP3.tif", "DP4/SoilMoisture_CECRSPPP_DP4.tif", "DP5/SoilMoisture_CECRSPPP_DP5.tif", "DP6/SoilMoisture_CECRSPPP_DP6.tif")
canp <- c("DP1/CanopyCover.tif", "DP2/CanopyCover_DP2.tif", "DP3/CanopyCover_CECRSPPP_DP3.tif", "DP4/CanopyCover_CECRSPPP_DP4.tif", "DP5/CanopyCover_CECRSPPP_DP5.tif", "DP6/CanopyCover_CECRSPPP_DP6.tif")
hydc <- c("DP1/Hydrochory.tif", "DP2/Hydrochory_DP2.tif", "DP3/Hydrochory_CECRSPPP_DP3.tif", "DP4/Hydrochory_CECRSPPP_DP4.tif", "DP5/Hydrochory_CECRSPPP_DP5.tif", "DP6/Hydrochory_CECRSPPP_DP6.tif")
prgs <- c("DP1/PropaguleSupply.tif", "DP2/PropaguleSupply_DP2.tif", "DP3/PropaguleSupply_CECRSPPP_DP3.tif", "DP4/PropaguleSupply_CECRSPPP_DP4.tif", "DP5/PropaguleSupply_CECRSPPP_DP5.tif", "DP6/PropaguleSupply_CECRSPPP_DP6.tif")
zooc <-c("DP1/Zoochory.tif", "DP2/Zoochory_DP2.tif", "DP3/Zoochory_CECRSPPP_DP3.tif", "DP4/Zoochory_CECRSPPP_DP4.tif", "DP5/Zoochory_CECRSPPP_DP5.tif", "DP6/Zoochory_CECRSPPP_DP6.tif")

## Define detection periods

detp <- 1:6

## Source file containing validation functions
## Note, the [...] should point to the sensitivity_analysis folder

source("[...]/sensitivity_functions.R")

## Loop through files and calculate Continuous Boyce Index
detections1<- readOGR("../Detections/DP1_200809.shp")
detections2 <- readOGR("../Detections/DP2_200910.shp")
detections3 <- readOGR("../Detections/DP3_201011.shp")
detections4 <- readOGR("../Detections/DP4_201112.shp")
detections5 <- readOGR("../Detections/DP5_201213.shp")
detections6 <- readOGR("../Detections/DP6_201314.shp")
detections7 <- readOGR("../Detections/DP7_201415.shp")
detections8 <- readOGR("../Detections/DP8_201516.shp")
detections9 <- readOGR("../Detections/DP9_201617.shp")
dets_1 <- rbind(
  detections2,
  detections3,
  detections4
)
dets_2 <- rbind(
  detections3,
  detections4,
  detections5
)
dets_3 <- rbind(
  detections4,
  detections5,
  detections6
)
dets_4 <- rbind(
  detections5,
  detections6,
  detections7
)
dets_5 <- rbind(
  detections6,
  detections7,
  detections8
)
dets_6 <- rbind(
  detections7,
  detections8,
  detections9
)
dets_ <- list(dets_1, dets_2, dets_3, dets_4, dets_5, dets_6)
for(i in 1:6){
  dets_[[i]]$V <- 1
}
cbis_ <- vector("list", 6)

## Computation of goodness of fit statistics

for(i in detp){
  # Loop through detection periods
  elev_i <- elev[[i]]
  rain_i <- rain[[i]]
  vegl_i <- vegl[[i]]
  soil_i <- soil[[i]]
  canp_i <- canp[[i]]
  hydc_i <- hydc[[i]]
  prgs_i <- prgs[[i]]
  zooc_i <- zooc[[i]]
  # Loop through resolutions
  cbis_[[i]] <- numeric(length(1:50))
  for(j in 1:50){
    cbis_[[i]][j] <- compute_bn_cbi(
      persistence = c(canp_i, soil_i, vegl_i),
      establishment = c(elev_i, rain_i),
      propagule = c(hydc_i, prgs_i, zooc_i), 
      detections = dets_[[i]],
      abundance_column = "V",
      persistence_wts = c(2, 1, 3), # numeric vec
      persistence_sd = j + 5, # numeric
      establishment_wts = c(1, 2), # numeric vec
      establishment_sd = j + 5, # numeric
      propagule_wts = c(3, 2, 3), # numeric vec
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
    cbis_[[4]],
    cbis_[[5]],
    cbis_[[6]]
  ),
  detp = paste("Detection period:", rep(detp, each = 50)),
  resl = 1:50 # input standard deviation
)

ggplot(data = raster_data, aes(x = resl, y = cbis)) +
  geom_point(size = 2) +
  geom_smooth(se = FALSE) +
  facet_wrap(~ detp, nrow = 2) +
  labs(x = "Standard deviation at child nodes", y = "Continuous Boyce Index (CBI)") +
  xlim(0, NA) +
  theme_bw()

save.image("Mexican_Bean_Tree_SD.Rdata")
