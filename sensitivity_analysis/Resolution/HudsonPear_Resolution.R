### Script:  SENSITIVITY ANALYSIS
### Species: CYLIROSE (HUDSON PEAR)
### Varying: RESOLUTION
### Author:  ALAN RYOU PEARSE

rm(list = ls())

## Working directory
## The [...] should point to the sensitivity_analysis folder

setwd("[...]/Data/HudsonPear/Riskfactors_susceptibility")

## Find raster files we want to process

files <- dir(pattern = ".tif$", recursive = TRUE)

## Rasters

temp <- files[grep("Temperature", files)]
rain <- files[grep("RainAnnual", files)]
vegl <- files[grep("VegetationLandUse", files)]
soil <- files[grep("SoilTexture", files)]
agoc <- files[grep("Agochory", files)]
hydc <- files[grep("Hydrochory", files)]
prgs <- files[grep("Propagule", files)]
zooc <- files[grep("Zoochory", files)]

## Split up by detection period

detp <- 1:4
cbis_ <- temp_ <- rain_ <- vegl_ <- soil_ <- agoc_ <- hydc_ <- prgs_ <- zooc_ <- vector("list", 4)
for(i in detp){
  temp_[[i]] <- temp[grep(paste0("DP", i), temp)]
  rain_[[i]] <- rain[grep(paste0("DP", i), rain)]
  vegl_[[i]] <- vegl[grep(paste0("DP", i), vegl)]
  soil_[[i]] <- soil[grep(paste0("DP", i), soil)]
  agoc_[[i]] <- agoc[grep(paste0("DP", i), agoc)]
  hydc_[[i]] <- hydc[grep(paste0("DP", i), hydc)]
  prgs_[[i]] <- prgs[grep(paste0("DP", i), prgs)]
  zooc_[[i]] <- zooc[grep(paste0("DP", i), zooc)]
}

## Source file containing validation functions
## Note, the [...] has to point to the sensitivity_analysis folder

source("[...]/sensitivity_functions.R")

## Loop through files and calculate Continuous Boyce Index
## This code assumes the working directory has already been set correctly and the folder structure from the GitHub repo has been preserved
indices <- sort(as.character(1:20))
ind_ord <- order(as.numeric(indices))
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

## Computation of goodness of fit statistics

for(i in detp){
  # Loop through detection periods
  temp_i <- temp_[[i]]
  rain_i <- rain_[[i]]
  vegl_i <- vegl_[[i]]
  soil_i <- soil_[[i]]
  agoc_i <- agoc_[[i]]
  hydc_i <- hydc_[[i]]
  prgs_i <- prgs_[[i]]
  zooc_i <- zooc_[[i]]
  # Loop through resolutions
  cbis_[[i]] <- numeric(20)
  cbis_[[i]][1] <- compute_bn_cbi(
    persistence = c(temp_i[1], rain_i[1]),
    establishment = c(vegl_i[1], soil_i[1]),
    propagule = c(agoc_i[1], hydc_i[1], prgs_i[1], zooc_i[1]), 
    detections = dets_[[i]],
    abundance_column = "V",
    persistence_wts = c(2, 2), # numeric vec
    persistence_sd = 15, # numeric
    establishment_wts = c(2, 1), # numeric vec
    establishment_sd = 15, # numeric
    propagule_wts = c(2, 3, 3, 2), # numeric vec
    propagule_sd = 15, # numeric
    suitability_sd = 10, # numeric
    susceptibility_sd = 10 # numeric
  )$Spearman.cor
  for(j in 2:20){
    cbis_[[i]][j] <- compute_bn_cbi(
      persistence = c(temp_i[j], rain_i[j]),
      establishment = c(vegl_i[j], soil_i[j]),
      propagule = c(agoc_i[j], hydc_i[j], prgs_i[j], zooc_i[j]), 
      detections = dets_[[i]],
      abundance_column = "V",
      persistence_wts = c(2, 2), # numeric vec
      persistence_sd = 15, # numeric
      establishment_wts = c(2, 1), # numeric vec
      establishment_sd = 15, # numeric
      propagule_wts = c(2, 3, 3, 2), # numeric vec
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
    cbis_[[4]]
  ),
  detp = paste("Detection period:", rep(detp, each = 20)),
  resl = (200 * 1:20) ^ 2 * 1e-4 # to hectares
)

ggplot(data = raster_data, aes(x = resl, y = cbis)) +
  geom_point(size = 2) +
  geom_smooth(se = FALSE) +
  facet_wrap(~ detp, nrow = 1) +
  labs(x = "Resolution (ha)", y = "Continuous Boyce Index (CBI)") +
  xlim(0, NA) +
  theme_bw()

save.image("HUDSON_PEAR_RESOLUTION.Rdata")