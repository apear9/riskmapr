expand_grid <- function(n_var, numerical_states = seq(0, 100, 25)){
  
  # Assign variable
  V <- numerical_states 
  
  # Grid expansion
  code_block <- "expand.grid("
  additions <- "V"
  if(n_var > 1){
    for(i in 1:(n_var - 1)){
      additions <- paste(additions, "V", sep = ", ")
    }
  }
  code_block <- paste0(code_block, additions, ")") # close bracket
  
  # Result
  return(eval(parse(text = code_block)))
  
}

cut_down_grid <- function(grid){
  
  # Input checking
  nc <- ncol(grid)
  nr <- nrow(grid)
  if(nc > 1){
    
    # Cut down grid
    ind <- rep(TRUE, nr)
    for(i in 2:nc){
      ind <- ind & (unlist(grid[, i-1]) < unlist(grid[, i]))
    }
    
    # Result
    return(grid[ind, ])
    
  } else {
    
    return(grid)
    
  }
  
}

replace <- function(r, grid){
  
  # Get unique values in raster in order
  uv <- sort(unique(r)) 
  r_copy <- r
  
  # Input checking
  nr <- nrow(grid)
  nc <- ncol(grid)
  if(nc != length(uv)) stop("The number of unique values in the raster r must match the number of of columns in grid.")
  
  # Replace raster values
  file_original <- filename(r)
  file_appendages <- paste0("_Grid_Levels_", 1:nr, ".tif")
  for(i in 1:nr){
    file_new <- gsub(".tif", file_appendages[i], file_original)
    for(j in 1:nc){
      r_copy[r[] == uv[j]] <- unlist(grid[i, j])
    }
    writeRaster(r_copy, file_new, overwrite = TRUE)
  }
  
  # When done, return nothing
  message("Processing complete.")
  invisible()
  
}