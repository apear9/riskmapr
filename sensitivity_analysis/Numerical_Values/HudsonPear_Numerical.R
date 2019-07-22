### Script:  SENSITIVITY ANALYSIS
### Species: CYLIROSE (MEXICAN BEAN TREE)
### Varying: NUMERICAL STATES
### Notes:   N/A
### Author:  ALAN RYOU PEARSE

rm(list = ls())

## Working directory
## The [...] should point to the sensitivity_analysis folder

setwd("[...]/Data/HudsonPear/Riskfactors_susceptibility")

## Find raster files we want to process

files <- dir(pattern = ".tif", recursive = TRUE)

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
## The [...] should point to sensitivity_analysis folder.

source("[...]/sensitivity_functions.R")

## Import detection records

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
## Set up validation for each risk factor 

# Elevation
cbis_temp <- vector("list", max(detp))
t1 <- Sys.time()
for(i in detp){
  n_files <- length(temp_[[i]])
  cbis_temp[[i]] <- numeric(n_files - 1)
  for(j in 2:n_files){
    cbis_temp[[i]][j-1] <- compute_bn_cbi(
      persistence = c(temp_[[i]][j], rain_[[i]][1]),
      establishment = c(vegl_[[i]][1], soil_[[i]][1]),
      propagule = c(agoc_[[i]][1], hydc_[[i]][1], prgs_[[i]][1], zooc_[[i]][1]), 
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
}
t2 <- Sys.time()
t2 - t1

## Rain 
cbis_rain <- vector("list", max(detp))
for(i in detp){
  n_files <- length(rain_[[i]])
  cbis_rain[[i]] <- numeric(n_files - 1)
  for(j in 2:n_files){
    cbis_rain[[i]][j-1] <- compute_bn_cbi(
      persistence = c(temp_[[i]][1], rain_[[i]][j]),
      establishment = c(vegl_[[i]][1], soil_[[i]][1]),
      propagule = c(agoc_[[i]][1], hydc_[[i]][1], prgs_[[i]][1], zooc_[[i]][1]), 
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
}

## Canopy 
cbis_agoc <- vector("list", max(detp))
for(i in detp){
  n_files <- length(agoc_[[i]])
  cbis_agoc[[i]] <- numeric(n_files - 1)
  for(j in 2:n_files){
    cbis_agoc[[i]][j-1] <- compute_bn_cbi(
      persistence = c(temp_[[i]][1], rain_[[i]][1]),
      establishment = c(vegl_[[i]][1], soil_[[i]][1]),
      propagule = c(agoc_[[i]][j], hydc_[[i]][1], prgs_[[i]][1], zooc_[[i]][1]), 
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
}

## soil 
cbis_soil <- vector("list", max(detp))
t1 <- Sys.time()
for(i in detp){
  n_files <- length(soil_[[i]])
  cbis_soil[[i]] <- numeric(n_files - 1)
  for(j in 2:n_files){
    cbis_soil[[i]][j-1] <- compute_bn_cbi(
      persistence = c(temp_[[i]][1], rain_[[i]][1]),
      establishment = c(vegl_[[i]][1], soil_[[i]][j]),
      propagule = c(agoc_[[i]][1], hydc_[[i]][1], prgs_[[i]][1], zooc_[[i]][1]), 
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
}

## vegl 
cbis_vegl <- vector("list", max(detp))
for(i in detp){
  n_files <- length(vegl_[[i]])
  cbis_vegl[[i]] <- numeric(n_files - 1)
  for(j in 2:n_files){
    cbis_vegl[[i]][j-1] <- compute_bn_cbi(
      persistence = c(temp_[[i]][1], rain_[[i]][1]),
      establishment = c(vegl_[[i]][j], soil_[[i]][1]),
      propagule = c(agoc_[[i]][1], hydc_[[i]][1], prgs_[[i]][1], zooc_[[i]][1]), 
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
}

## hydc 
cbis_hydc <- vector("list", max(detp))
for(i in detp){
  n_files <- length(hydc_[[i]])
  cbis_hydc[[i]] <- numeric(n_files - 1)
  for(j in 2:n_files){
    cbis_hydc[[i]][j-1] <- compute_bn_cbi(
      persistence = c(temp_[[i]][1], rain_[[i]][1]),
      establishment = c(vegl_[[i]][1], soil_[[i]][1]),
      propagule = c(agoc_[[i]][1], hydc_[[i]][j], prgs_[[i]][1], zooc_[[i]][1]), 
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
}

## prgs 
cbis_prgs <- vector("list", max(detp))
for(i in detp){
  n_files <- length(prgs_[[i]])
  cbis_prgs[[i]] <- numeric(n_files - 1)
  for(j in 2:n_files){
    cbis_prgs[[i]][j-1] <- compute_bn_cbi(
      persistence = c(temp_[[i]][1], rain_[[i]][1]),
      establishment = c(vegl_[[i]][1], soil_[[i]][1]),
      propagule = c(agoc_[[i]][1], hydc_[[i]][1], prgs_[[i]][j], zooc_[[i]][1]), 
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
}

## zooc 
cbis_zooc <- vector("list", max(detp))
for(i in detp){
  n_files <- length(zooc_[[i]])
  cbis_zooc[[i]] <- numeric(n_files - 1)
  for(j in 2:n_files){
    cbis_zooc[[i]][j-1] <- compute_bn_cbi(
      persistence = c(temp_[[i]][1], rain_[[i]][1]),
      establishment = c(vegl_[[i]][1], soil_[[i]][1]),
      propagule = c(agoc_[[i]][1], hydc_[[i]][1], prgs_[[i]][1], zooc_[[i]][j]), 
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
}

# Save data
save.image("HUDSON_PEAR_NUMERICAL.Rdata")