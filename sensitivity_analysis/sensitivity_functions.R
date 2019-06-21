### Script: SENSITIVITY ANALYSIS FUNCTIONS
### Author: ALAN RYOU PEARSE

### Sensitivity analysis -- parameters, inputs, and data resolution
### Alan R. Pearse, 29/04/2019
### With Jens G. Froese and Grant S. Hamilton

# Load packages
library(rgdal)
library(shiny)
library(raster)
library(sp)
library(truncnorm)
library(stringr)

# Define functions for computing the moments of distributions
exp_discrete <- function(x){
  p <- x
  s <- as.numeric(names(x))
  sum(p * s)
}
ex2_discrete <- function(x){
  p <- x
  s <- as.numeric(names(x))
  sum(p * (s^2))
}
std_discrete <- function(x){
  sqrt(ex2_discrete(x) - exp_discrete(x)^2)
}
# Inputs:
# persistence, persistence_weights, persistence_sd, 
# establishment, establishment_weights, establishment_sd,
# propagule, propagule_weights, propagule_sd,
# suitability_sd, susceptibility_sd
get_predictions <- function(
  persistence, # list
  persistence_wts, # numeric vec 
  persistence_sd, # numeric
  establishment, # list
  establishment_wts, # numeric vec
  establishment_sd, # numeric 
  propagule, # list
  propagule_wts, # numeric vec
  propagule_sd, # numeric
  suitability_sd, # numeric
  susceptibility_sd # numeric
){
  # Subsets as required for the analysis
  per_wets <- persistence_wts 
  est_wets <- establishment_wts 
  prg_wets <- propagule_wts
  # Empty numeric vectors
  st <- st_sd <- sc <- sc_sd <- numeric(length(persistence))
  # Main loop
  for(i in 1:length(persistence)){
    # Establishment
    est_vars <- establishment[[i]]
    est_mean <- sum(est_vars * est_wets)/sum(est_wets)
    est <- dtruncnorm(seq(0, 100, 25), 0, 100, est_mean, establishment_sd)
    est <- est/sum(est)
    names(est) <- seq(0, 100, 25)
    # Persistence
    per_vars <- persistence[[i]]
    per_mean <- sum(per_vars * per_wets)/sum(per_wets)
    per <- dtruncnorm(seq(0, 100, 25), 0, 100, per_mean, persistence_sd)
    per <- per/sum(per)
    names(per) <- seq(0, 100, 25)
    # Suitability by marginalisation using law of total probability
    n_est <- length(est)
    n_per <- length(per)
    j_mat <- matrix(0, nrow = n_est * n_per, ncol = 5)
    cnt <- 1
    est_x <- as.numeric(names(est))
    per_x <- as.numeric(names(per))
    for(j in 1:n_est){
      for(k in 1:n_per){
        p_jk <- dtruncnorm(seq(0, 100, 25), 0, 100, (sum(est_wets) * est_x[j] + sum(per_wets) * per_x[k])/(sum(est_wets) + sum(per_wets)), suitability_sd)
        j_mat[cnt, ] <- p_jk/sum(p_jk) * est[j] * per[k]
        cnt <- cnt + 1
      }
    }
    suit <- colSums(j_mat)
    names(suit) <- seq(0, 100, 25)
    suit_wets <- sum(est_wets) + sum(per_wets)
    # Take expectation as prediction
    st[i] <- exp_discrete(suit)
    st_sd[i] <- std_discrete(suit)
    # Now construct the propagule pressure node
    prg_vars <- propagule[[i]]
    prg_mean <- sum(prg_vars * prg_wets)/sum(prg_wets)
    prg <- dtruncnorm(seq(0, 100, 25), 0, 100, prg_mean, propagule_sd)
    prg <- prg/sum(prg)
    names(prg) <- seq(0, 100, 25)
    # Now construct the susceptibility node
    n_prg <- length(prg)
    n_sut <- length(suit)
    j_mat <- matrix(0, nrow = n_sut * n_prg, ncol = 5)
    cnt <- 1
    suit_x <- as.numeric(names(suit))
    prg_x <- as.numeric(names(prg))
    for(j in 1:n_sut){
      for(k in 1:n_prg){
        p_jk <- dtruncnorm(seq(0, 100, 25), 0, 100, (sum(suit_wets) * suit_x[j] + sum(prg_wets) * prg_x[k])/(sum(suit_wets) + sum(prg_wets)), susceptibility_sd)
        j_mat[cnt, ] <- p_jk/sum(p_jk) * suit[j] * prg[k]
        cnt <- cnt + 1
      }
    }
    susc <- colSums(j_mat)
    names(susc) <- seq(0, 100, 25)
    # Take expectation as the prediction
    sc[i] <- exp_discrete(susc)
    sc_sd[i] <- std_discrete(susc)
    
  }
  # Return
  list(
    Suitability = st, Suitability_SD = st_sd, 
    Susceptibility = sc, Susceptibility_SD = sc_sd
  )
  
}
get_data <- function(
  persistence,
  establishment,
  propagule, 
  return_as = "list",
  return_raster = TRUE
){
  # Create column indices
  i_per <- 1:length(persistence)
  i_est <- 1:length(establishment) + length(i_per)
  i_prg <- 1:length(propagule) + length(i_per) + length(i_est)
  # Read in rasters as a stack
  raster_data <- stack(c(persistence, establishment, propagule))
  # Get in data frame format
  raster_data_df <- dplyr::distinct(
    as.data.frame(raster_data)
  )
  # Remove NAs
  raster_data_df <- na.omit(raster_data_df)
  ind <- apply(raster_data_df, 1, function(x) any(x < 0 | x > 100))
  raster_data_df <- raster_data_df[!ind, ]
  # Return in format requested by user (data.frame or as list)
  if(return_as == "data.frame"){
    return(raster_data_df)
  } else {
    n <- nrow(raster_data_df)
    per_list <- est_list <- prg_list <- vector("list", n)
    for(i in 1:n){
      per_list[[i]] <- raster_data_df[i, i_per]
      est_list[[i]] <- raster_data_df[i, i_est]
      prg_list[[i]] <- raster_data_df[i, i_prg]
    }
    to_return <- list(
      persistence = per_list,
      establishment = est_list,
      propagule = prg_list
    )
    if(return_raster) to_return$rasters <- raster_data
    return(
      to_return
    )
  }
}
# Function to construct background list for Continuous Boyce Index calculations
construct_raster_list <- function(raster_layer){
  raster_table <- table(raster_layer[])
  raster_list  <- rep(
    as.numeric(names(raster_table)),
    raster_table
  )
  return(sort(raster_list))
}
# Function to construct observed list for Continuous Boyce Index calculations
construct_detections_list <- function(raster_layer, detections, abundance_column){
  SI <- extract(raster_layer, detections)
  SI[is.na(SI)] <- 10 # assign low value to detections outside of risk area
  abundances <- unlist(detections@data[, abundance_column])
  vals <- unique(SI)
  n <- length(vals)
  nums <- numeric(n)
  for(i in 1:n){
    ind <- SI == vals[i]
    abd <- sum(abundances[ind])
    nums[i] <- floor(abd)
  }
  SIs <- rep(vals, nums)
  list(sort(SIs))
}
# Function to compute the Continuous Boyce Index for predictions from the BBN
compute_bn_cbi <- function(
  persistence,
  establishment,
  propagule,
  detections,
  abundance_column = "Abundance",
  ... # Additional arguments to either get_predictions
){
  # Extract out data from the rasters
  dat <- get_data(persistence, establishment, propagule)
  # Get a single RasterLayer to store predictions
  raster_layer <- dat$rasters[[1]]
  # Get predictions 
  prd <- get_predictions(
    persistence = dat$persistence,
    establishment = dat$establishment, 
    propagule = dat$propagule, 
    ...
  )
  # Get raster as data.frame
  df_all <- as.data.frame(
    dat$rasters
  )
  names(df_all) <- paste0("V", 1:ncol(df_all))
  
  # Get unique combinations
  df_unique <- as.data.frame(
    cbind(
      do.call(rbind, dat$persistence),
      do.call(rbind, dat$establishment),
      do.call(rbind, dat$propagule)
    )
  )
  names(df_unique) <- paste0("V", 1:ncol(df_unique))
  df_unique$Vals <- prd$Susceptibility
  # Perform a join
  df_all <- dplyr::left_join(df_all, df_unique, by = paste0("V", 1:ncol(df_all)))
  # Store values in raster
  raster_layer[] <- floor(df_all$Vals)
  # Set parameters for CBI computation
  nclass = 0 # defaults to moving window (continuous, classification-independent) computation with arguments
  res = 100 # resolution factor (i.e. 100 computations across model-predicted value range)
  window.w = 15 # moving window width
  PEplot = FALSE # no PEplot is generated 
  ras_list <- construct_raster_list(raster_layer)
  det_list <- construct_detections_list(raster_layer, detections, abundance_column)[[1]]
  # Compute CBI
  cbi <- ecospat:::ecospat.boyce(
    ras_list,
    det_list,
    nclass, 
    window.w, 
    res, 
    PEplot
  )
  # Return
  return(cbi)
}
