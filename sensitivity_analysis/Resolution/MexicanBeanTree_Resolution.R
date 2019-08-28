### Script:  SENSITIVITY ANALYSIS
### Species: CECROPIA (MEXICAN BEAN TREE)
### Varying: RESOLUTION
### Notes:   N/A
### Author:  ALAN RYOU PEARSE

rm(list = ls())

## Working directory
## The [...] should point to the folder sensitivity_analysis

setwd("[...]/Data/MexicanBeanTree/Riskfactors_susceptibility")

## Find raster files we want to process

files <- dir(pattern = ".tif$", recursive = TRUE)

## Rasters

elev <- files[grep("Elevation", files)]
canp <- files[grep("CanopyCover", files)]
rain <- files[grep("RainDriestMonth", files)]
vegl <- files[grep("VegetationLandUse", files)]
soil <- files[grep("SoilMoisture", files)]
hydc <- files[grep("Hydrochory", files)]
prgs <- files[grep("Propagule", files)]
zooc <- files[grep("Zoochory", files)]

## Split up by detection period

detp <- 1:6
cbis_ <- elev_ <- rain_ <- vegl_ <- soil_ <- canp_ <- hydc_ <- prgs_ <- zooc_ <- vector("list", 6)
for(i in detp){
  elev_[[i]] <- elev[grep(paste0("DP", i), elev)]
  rain_[[i]] <- rain[grep(paste0("DP", i), rain)]
  vegl_[[i]] <- vegl[grep(paste0("DP", i), vegl)]
  soil_[[i]] <- soil[grep(paste0("DP", i), soil)]
  canp_[[i]] <- canp[grep(paste0("DP", i), canp)]
  hydc_[[i]] <- hydc[grep(paste0("DP", i), hydc)]
  prgs_[[i]] <- prgs[grep(paste0("DP", i), prgs)]
  zooc_[[i]] <- zooc[grep(paste0("DP", i), zooc)]
}

## Source file containing validation functions
## Note, the [...] has to point to the sensitivity_analysis folder

source("[...]/sensitivity_functions.R")

## Loop through files and calculate Continuous Boyce Index
indices <- sort(as.character(1:40))
ind_ord <- order(as.numeric(indices))

## Note, this code assumes the working directory has been correctly set
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

## Computation of goodness of fit statistics

for(i in detp){
  # Loop through detection periods
  elev_i <- elev_[[i]]
  rain_i <- rain_[[i]]
  vegl_i <- vegl_[[i]]
  soil_i <- soil_[[i]]
  canp_i <- canp_[[i]]
  hydc_i <- hydc_[[i]]
  prgs_i <- prgs_[[i]]
  zooc_i <- zooc_[[i]]
  # Loop through resolutions
  cbis_[[i]] <- numeric(40)
  cbis_[[i]][1] <- compute_bn_cbi(
    persistence = c(canp_i[1], soil_i[1], vegl_i[1]),
    establishment = c(elev_i[1], rain_i[1]),
    propagule = c(hydc_i[1], prgs_i[1], zooc_i[1]), 
    detections = dets_[[i]],
    abundance_column = "V",
    persistence_wts = c(2, 1, 3), # numeric vec
    persistence_sd = 15, # numeric
    establishment_wts = c(1, 2), # numeric vec
    establishment_sd = 15, # numeric
    propagule_wts = c(2, 3, 3), # numeric vec
    propagule_sd = 15, # numeric
    suitability_sd = 10, # numeric
    susceptibility_sd = 10 # numeric
  )$Spearman.cor
  for(j in 2:40){
    cbis_[[i]][j] <- compute_bn_cbi(
      persistence = c(canp_i[j], soil_i[j], vegl_i[j]),
      establishment = c(elev_i[j], rain_i[j]),
      propagule = c(hydc_i[j], prgs_i[j], zooc_i[j]), 
      detections = dets_[[i]],
      abundance_column = "V",
      persistence_wts = c(2, 1, 3), # numeric vec
      persistence_sd = 15, # numeric
      establishment_wts = c(1, 2), # numeric vec
      establishment_sd = 15, # numeric
      propagule_wts = c(2, 3, 3), # numeric vec
      propagule_sd = 15, # numeric
      suitability_sd = 10, # numeric
      susceptibility_sd = 10 # numeric
    )$Spearman.cor
  }
  cbis_[[i]] <- cbis_[[i]][ind_ord]
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
  detp = paste("Detection period:", rep(detp, each = 40)),
  resl = (100 * 1:40) ^ 2 * 1e-4 # to hectares
)

ggplot(data = raster_data, aes(x = resl, y = cbis)) +
  geom_point(size = 2) +
  geom_smooth(se = FALSE) +
  facet_wrap(~ detp, nrow = 2) +
  labs(x = "Resolution (ha)", y = "Continuous Boyce Index (CBI)") +
  xlim(0, NA) +
  theme_bw()

save.image("MEXICAN_BEAN_TREE_RESOLUTION.Rdata")
