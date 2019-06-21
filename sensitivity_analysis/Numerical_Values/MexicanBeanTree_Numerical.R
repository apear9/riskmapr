### Script:  SENSITIVITY ANALYSIS
### Species: CECROPIA (MEXICAN BEAN TREE)
### Varying: NUMERICAL STATES
### Notes:   N/A
### Author:  ALAN RYOU PEARSE

rm(list = ls())

## Working directory
## [...] should point to S4_Appendix 

setwd("[...]/MexicanBeanTree")

## Find raster files we want to process

files <- dir(pattern = ".tif", recursive = TRUE)

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
## [...] should point to sensitivity_analysis

source("[...]/sensitivity_functions.R")

## Import detection records
## [...] should point to S4_Appendix 

detections1<- readOGR("[...]/MexicanBeanTree/Detections/DP1_200809.shp")
detections2 <- readOGR("[...]/MexicanBeanTree/Detections/DP2_200910.shp")
detections3 <- readOGR("[...]/MexicanBeanTree/Detections/DP3_201011.shp")
detections4 <- readOGR("[...]/MexicanBeanTree/Detections/DP4_201112.shp")
detections5 <- readOGR("[...]/MexicanBeanTree/Detections/DP5_201213.shp")
detections6 <- readOGR("[...]/MexicanBeanTree/Detections/DP6_201314.shp")
detections7 <- readOGR("[...]/MexicanBeanTree/Detections/DP7_201415.shp")
detections8 <- readOGR("[...]/MexicanBeanTree/Detections/DP8_201516.shp")
detections9 <- readOGR("[...]/MexicanBeanTree/Detections/DP9_201617.shp")
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

## Set up validation for each risk factor 

# Elevation
cbis_elev <- vector("list", max(detp))
t1 <- Sys.time()
for(i in detp){
  n_files <- length(elev_[[i]])
  cbis_elev[[i]] <- numeric(n_files - 1)
  for(j in 2:n_files){
    cbis_elev[[i]][j-1] <- compute_bn_cbi(
      persistence = c(canp_[[i]][1], soil_[[i]][1], vegl_[[i]][1]),
      establishment = c(elev_[[i]][j], rain_[[i]][1]),
      propagule = c(hydc_[[i]][1], prgs_[[i]][1], zooc_[[i]][1]), 
      detections = dets_[[i]],
      abundance_column = "V",
      persistence_wts = c(2, 1, 3), # numeric vec
      persistence_sd = 15, # numeric
      establishment_wts = c(1, 2), # numeric vec
      establishment_sd = 15, # numeric
      propagule_wts = c(3, 2, 3), # numeric vec
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
      persistence = c(canp_[[i]][1], soil_[[i]][1], vegl_[[i]][1]),
      establishment = c(elev_[[i]][1], rain_[[i]][j]),
      propagule = c(hydc_[[i]][1], prgs_[[i]][1], zooc_[[i]][1]), 
      detections = dets_[[i]],
      abundance_column = "V",
      persistence_wts = c(2, 1, 3), # numeric vec
      persistence_sd = 15, # numeric
      establishment_wts = c(1, 2), # numeric vec
      establishment_sd = 15, # numeric
      propagule_wts = c(3, 2, 3), # numeric vec
      propagule_sd = 15, # numeric
      suitability_sd = 10, # numeric
      susceptibility_sd = 10 # numeric
    )$Spearman.cor
  }
}

## Canopy 
cbis_canp <- vector("list", max(detp))
for(i in detp){
  n_files <- length(canp_[[i]])
  cbis_canp[[i]] <- numeric(n_files - 1)
  for(j in 2:n_files){
    cbis_canp[[i]][j-1] <- compute_bn_cbi(
      persistence = c(canp_[[i]][j], soil_[[i]][1], vegl_[[i]][1]),
      establishment = c(elev_[[i]][1], rain_[[i]][1]),
      propagule = c(hydc_[[i]][1], prgs_[[i]][1], zooc_[[i]][1]), 
      detections = dets_[[i]],
      abundance_column = "V",
      persistence_wts = c(2, 1, 3), # numeric vec
      persistence_sd = 15, # numeric
      establishment_wts = c(1, 2), # numeric vec
      establishment_sd = 15, # numeric
      propagule_wts = c(3, 2, 3), # numeric vec
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
      persistence = c(canp_[[i]][1], soil_[[i]][j], vegl_[[i]][1]),
      establishment = c(elev_[[i]][1], rain_[[i]][1]),
      propagule = c(hydc_[[i]][1], prgs_[[i]][1], zooc_[[i]][1]), 
      detections = dets_[[i]],
      abundance_column = "V",
      persistence_wts = c(2, 1, 3), # numeric vec
      persistence_sd = 15, # numeric
      establishment_wts = c(1, 2), # numeric vec
      establishment_sd = 15, # numeric
      propagule_wts = c(3, 2, 3), # numeric vec
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
      persistence = c(canp_[[i]][1], soil_[[i]][1], vegl_[[i]][j]),
      establishment = c(elev_[[i]][1], rain_[[i]][1]),
      propagule = c(hydc_[[i]][1], prgs_[[i]][1], zooc_[[i]][1]), 
      detections = dets_[[i]],
      abundance_column = "V",
      persistence_wts = c(2, 1, 3), # numeric vec
      persistence_sd = 15, # numeric
      establishment_wts = c(1, 2), # numeric vec
      establishment_sd = 15, # numeric
      propagule_wts = c(3, 2, 3), # numeric vec
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
      persistence = c(canp_[[i]][1], soil_[[i]][1], vegl_[[i]][1]),
      establishment = c(elev_[[i]][1], rain_[[i]][1]),
      propagule = c(hydc_[[i]][j], prgs_[[i]][1], zooc_[[i]][1]), 
      detections = dets_[[i]],
      abundance_column = "V",
      persistence_wts = c(2, 1, 3), # numeric vec
      persistence_sd = 15, # numeric
      establishment_wts = c(1, 2), # numeric vec
      establishment_sd = 15, # numeric
      propagule_wts = c(3, 2, 3), # numeric vec
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
      persistence = c(canp_[[i]][1], soil_[[i]][1], vegl_[[i]][1]),
      establishment = c(elev_[[i]][1], rain_[[i]][1]),
      propagule = c(hydc_[[i]][1], prgs_[[i]][j], zooc_[[i]][1]), 
      detections = dets_[[i]],
      abundance_column = "V",
      persistence_wts = c(2, 1, 3), # numeric vec
      persistence_sd = 15, # numeric
      establishment_wts = c(1, 2), # numeric vec
      establishment_sd = 15, # numeric
      propagule_wts = c(3, 2, 3), # numeric vec
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
      persistence = c(canp_[[i]][1], soil_[[i]][1], vegl_[[i]][1]),
      establishment = c(elev_[[i]][1], rain_[[i]][1]),
      propagule = c(hydc_[[i]][1], prgs_[[i]][1], zooc_[[i]][j]), 
      detections = dets_[[i]],
      abundance_column = "V",
      persistence_wts = c(2, 1, 3), # numeric vec
      persistence_sd = 15, # numeric
      establishment_wts = c(1, 2), # numeric vec
      establishment_sd = 15, # numeric
      propagule_wts = c(3, 2, 3), # numeric vec
      propagule_sd = 15, # numeric
      suitability_sd = 10, # numeric
      susceptibility_sd = 10 # numeric
    )$Spearman.cor
  }
}

# Save data
save.image("MEXICAN_BEAN_TREE_NUMERICAL.Rdata")