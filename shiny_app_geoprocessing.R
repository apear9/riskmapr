library(rgdal)
library(raster)
library(rgeos)
library(sp)
library(zip)

ui <- fluidPage(
  
  # App title
  
  titlePanel("Tools for converting weed detection records into spatial risk proxies"),
  
  # Sidebar panel for inputs ----
  
  sidebarLayout(
    
    sidebarPanel(
      
      selectInput(
        inputId = "which",
        label = "Select the geoprocessing workflow you need:",
        choices = c(
          "None", 
          "Project shapefile", 
          #"Project raster", 
          "Crop to extent", 
          "Propagule supply", 
          "Anemochory", 
          "Zoochory", 
          "Agochory", 
          "Hydrochory", 
          "Stream edge to raster"#, 
          #"Recode stream raster"
          ), 
        selected = "None", 
        multiple = FALSE
      ),
      
      uiOutput("Stream_raster"),
      
      uiOutput("Stream_value"),
      
      # uiOutput("Other_value"),
      
      uiOutput("Stream_shapefile"),
      
      uiOutput("Generic_shapefile"),
      
      uiOutput("Generic_raster"),
      
      uiOutput("Detections"),
      
      uiOutput("Reference_raster"),
      
      uiOutput("Column"),
      
      uiOutput("Abundance_thresholds"),
      
      uiOutput("Distance_thresholds"),
      
      uiOutput("Proxy_levels"),
      
      uiOutput("Max_radius"),
      
      uiOutput("Proj4string"),
      
      uiOutput("Output_name"),
      
      uiOutput("submit_button"),
      
      uiOutput("download_button")
      
    ),
    
    mainPanel(
      plotOutput("plot")
    )
  )
  
)

server <- function(input, output){
  
  # Function to convert stream edges to a raster format
  stream_edges_to_raster <- function(stream_edges, reference_raster){
    
    reference_raster[] <- 1
    mask(reference_raster, stream_edges)
    
  }
  
  # Function to create hydrochory raster
  hydrochory <- function(stream_raster, detections, reference_raster, distance_thresholds, proxy_levels, max_radius = max(distance_thresholds)){
    
    # Check whether the number of distance thresholds and proxy levels agree
    n_p_lvl <- length(proxy_levels)
    n_d_tsd <- length(distance_thresholds)
    if(n_p_lvl + 1 != n_d_tsd){
      stop("The number of proxy levels and distance thresholds must be the same.")
    }
    
    # mask reference raster within maximum distance threshold
    mx_d_th <- max(max(distance_thresholds), max_radius)
    mask_rs <- mask_within_buffer(detections, reference_raster, max_radius)
    # stream_raster <- mask_within_buffer_streams(detections, stream_raster, max_radius)
    
    # Convert stream raster to SpatialPixelsDataFrame and extract coordinates
    stream_raster <- crop(stream_raster, extent(mask_rs))
    strm_px <- as(stream_raster, "SpatialPixelsDataFrame")
    strm_cd <- coordinates(strm_px)
    
    # Extract coordinates from detections
    dett_cd <- coordinates(detections)
    
    # Calculate distance from stream cells to detections
    n_coord <- nrow(strm_cd)
    n_detec <- nrow(detections@data)
    matrix_ <- matrix(0, nrow = n_coord, ncol = n_detec)
    for(i in 1:n_coord){
      matrix_[i,] <- euc_dist_vec(strm_cd[i,], dett_cd)
    }
    closest <- apply(matrix_, 1, min)
    
    # Compare distances to thresholds
    matrix_ <- matrix(0, nrow = n_coord, ncol = n_p_lvl)
    for(i in 1:n_p_lvl){
      matrix_[,i] <- proxy_levels[i] * 
        (is_between(closest, distance_thresholds[i], distance_thresholds[i+1]) * 1)
    }
    proxies <- apply(matrix_, 1, max)
    
    # Join to SpatialPixelsDataFrame
    strm_px@data$proxy <- proxies
    strm_px@data[] <- strm_px@data[, "proxy"]
    
    # Turn into raster
    hydrochory <- raster(strm_px)
    hydrochory[is.na(hydrochory[])] <- 0
    # hydrochory <- mask_within_buffer_streams(detections, hydrochory, max_radius)
    
    # Return
    return(hydrochory)
    
  }
  
  # Function to create propagule supply raster
  propagule_supply <- function(detections, reference_raster, column, abundance_thresholds, proxy_levels, max_radius){
    
    # Check whether the number of distance thresholds and proxy levels agree
    n_p_lvl <- length(proxy_levels)
    n_a_tsd <- length(abundance_thresholds)
    if(n_p_lvl + 1 != n_a_tsd){
      stop("The number of proxy levels and distance thresholds must be the same.")
    }
    
    # mask reference raster within maximum distance threshold
    mask_rs <- mask_within_buffer(detections, reference_raster, max_radius)
    
    # Turn this into a SpatialPixelsDataFrame
    mask_px <- as(mask_rs, "SpatialPixelsDataFrame")
    mask_cd <- coordinates(mask_px)
    dett_cd <- coordinates(detections)
    n_coord <- nrow(mask_cd)
    n_detec <- nrow(dett_cd)
    
    # Calculate some distances, calculate weed abundance
    matrix_ <- matrix(0, nrow = n_coord, ncol = n_detec)
    for(i in 1:n_coord){
      matrix_[i,] <- euc_dist_vec(mask_cd[i,], dett_cd)
    }
    n_detns <- num_detections_within(matrix_, detections, column, max_radius)
    
    # Match these numbers to the abundance thresholds for the proxies
    matrix_ <- matrix(0, n_coord, n_p_lvl)
    for(i in 1:n_p_lvl){
      matrix_[,i] <- proxy_levels[i] * 
        (is_between(n_detns, abundance_thresholds[i], abundance_thresholds[i+1]) * 1)
    }
    proxies <- apply(matrix_, 1, max)
    
    # Assign the proxies to the SpatialPixelsDataFrame
    mask_px@data$proxy <- proxies
    mask_px@data[] <- mask_px@data$proxy
    
    # Convert this to a raster layer
    propagule_supply <- raster(mask_px)
    
    # Return
    return(propagule_supply)
    
  }
  
  # Function to mask within buffer around detections
  mask_within_buffer <- function(detections, reference_raster, max_radius){
    
    # Check that all arguments are reasonable
    if(max_radius <= 0){
      stop("The argument max_radius must be greater than 0.")
    }
    
    # Buffer and mask
    buffered_detections <- gBuffer(detections, FALSE, width = max_radius)
    reference_raster <- crop(reference_raster, extent(buffered_detections))
    reference_raster[] <- 1
    masked_raster <- mask(reference_raster, mask = buffered_detections)
    return(masked_raster)
    
  }
  
  # Function to do the same as above but for stream rasters only
  mask_within_buffer_streams <- function(detections, reference_raster, max_radius){
    
    # Check that all arguments are reasonable
    if(max_radius <= 0){
      stop("The argument max_radius must be greater than 0.")
    }
    
    # Buffer and mask
    buffered_detections <- gBuffer(detections, FALSE, width = max_radius)
    reference_raster <- crop(reference_raster, extent(buffered_detections))
    masked_raster <- mask(reference_raster, mask = buffered_detections)
    return(masked_raster)
    
  }
  
  # Function to create a zoochory raster
  zoochory <- function(detections, reference_raster, distance_thresholds, proxy_levels, max_radius = max(distance_thresholds)){
    
    # Check whether the number of distance thresholds and proxy levels agree
    n_p_lvl <- length(proxy_levels)
    n_d_tsd <- length(distance_thresholds)
    if(n_p_lvl + 1 != n_d_tsd){
      stop("The number of proxy levels and distance thresholds must be the same.")
    }
    
    # mask reference raster within maximum distance threshold
    mx_d_th <- max(max(distance_thresholds), max_radius)
    mask_rs <- mask_within_buffer(detections, reference_raster, mx_d_th)
    
    # Turn this into a SpatialPixelsDataFrame
    mask_px <- as(mask_rs, "SpatialPixelsDataFrame")
    mask_cd <- coordinates(mask_px)
    dett_cd <- coordinates(detections)
    n_coord <- nrow(mask_cd)
    n_detec <- nrow(dett_cd)
    
    # Calculate some distances, apply thresholds, and take the maximum within the rows of the data.
    matrix_ <- matrix(0, nrow = n_coord, ncol = n_detec)
    for(i in 1:n_coord){
      matrix_[i,] <- euc_dist_vec(mask_cd[i,], dett_cd)
    }
    
    # Obtain minimum distance from source
    closest <- apply(matrix_, 1, min)
    
    # Check against distance thresholds
    matrix_ <- matrix(0, n_coord, n_p_lvl)
    for(i in 1:n_p_lvl){
      matrix_[,i] <- proxy_levels[i] * 
        (is_between(closest, distance_thresholds[i], distance_thresholds[i+1]) * 1)
    }
    proxies <- apply(matrix_, 1, max)
    
    # Join back to SpatialPixelsDataFrame
    mask_px@data$proxy <- proxies
    mask_px@data[] <- mask_px$proxy
    
    # Turn into raster
    zoochory <- raster(mask_px)
    
    # Return
    return(zoochory)
    
  }
  
  # Function to count the detections within a fixed radius
  num_detections_within <- function(matrix_, detections, column, distance_threshold){
    
    # Check that inputs are correct
    n_col <- length(column)
    if(n_col != 1){
      stop("Only one column can be specified for the argument column.")
    }
    
    # Calculate the number of detections within the region
    ((matrix_ <= distance_threshold) * 1) %*% detections@data[, column]
    
  }
  
  # Same as above but requiring no column argument
  num_detections_within_no_column <- function(matrix_, detections, distance_threshold){
    
    # Calculate the number of detections within the region
    sum((matrix_ <= distance_threshold) * 1)
    
  }
  
  # Convert stream rasters with two values (one for streams, another for others) to a raster with ones for stream cells and NAs elsewhere
  binary_stream_to_unary_stream <- function(stream_raster, stream_value){
    
    # Turn stream cells to 1
    if(stream_value != 1){
      stream_raster[stream_raster[] == stream_value] <- 1
    }
    
    # Turn non-stream cells to NA
    stream_raster[stream_raster[] != stream_value] <- NA
    
    # Return streams
    stream_raster
    
  }
  
  # Function to map thresholds to proxy levels
  is_between <- function(v, b1, b2){
    
    # Sort the bounds
    lower <- min(c(b1, b2))
    upper <- max(c(b1, b2))
    
    # Return logical for comparison
    lower <= v & v <= upper
    
  }
  
  # Function to measure Euclidean distance
  euc_dist <- function(x, y){
    
    d <- x - y
    sqrt(sum(d^2))
    
  }
  
  # Function to vectorise the above process
  euc_dist_vec <- function(x, y){
    
    apply(y, 1, euc_dist, x)
    
  }
  
  # This code updates the UI
  output$Stream_raster <- renderUI(
    {
      if(input$which %in% c("Hydrochory", "Agochory", "Recode stream raster")){
        fileInput("stream_raster", "Upload a raster of a linear feature", FALSE, ".tif")
      }
    }
  )
  
  output$Stream_value <- renderUI(
    {
      if(input$which == "Recode stream raster"){
        numericInput("stream_value", "Enter the value for the stream cells", 1, min = -1e9, max = 1e9)
      }
    }
  )
  
  # output$Other_value <- renderUI(
  #   {
  #     if(input$which == "Recode stream raster"){
  #       numericInput("other_value", "Enter the value for the non-stream cells", 0, min = -1e9, max = 1e9)
  #     }
  #   }
  # )
  
  output$Stream_shapefile <- renderUI(
    {
      if(input$which == "Stream edge to raster"){
        fileInput("stream_shapefile", "Upload a shapefile of streams", TRUE)
      }
    }
  )
  
  output$Detections <- renderUI(
    {
      if(input$which %in% c("Propagule supply", "Zoochory", "Anemochory", "Hydrochory", "Agochory", "Crop to extent")){
        fileInput("detections", "Upload shapefile of the detection records", TRUE)
      }
    }
  )
  
  output$Generic_raster <- renderUI(
    {
      if(input$which %in% c("Crop to extent", "Project raster", "Project shapefile")){
        fileInput("generic_raster", "Upload a raster (.tif extension)", FALSE, ".tif")
      }
    }
  )
  
  output$Generic_shapefile <- renderUI(
    {
      if(input$which %in% c("Project shapefile")){
        fileInput("generic_shapefile", "Upload a shapefile", TRUE)
      }
    }
  )
  
  # output$Proj4string <- renderUI(
  #   {
  #     if(input$which %in% c("Project shapefile", "Project raster")){
  #       textInput("proj4string", "Enter a proj4string (see spatialreference.org for proj4strings):", "+init=epsg:4326")
  #     }
  #   }
  # )
  
  output$Reference_raster <- renderUI(
    {
      if(input$which %in% c("Propagule supply", "Zoochory", "Anemochory", "Hydrochory", "Zoochory", "Agochory", "Stream edge to raster")){
        fileInput("reference_raster", "Upload the reference raster", FALSE, ".tif")
      }
    }
  )
  
  output$Column <- renderUI(
    {
      if(input$which == "Propagule supply"){
        textInput("column", "Name of the column containing the number of detections per record", "")
      }
    }
  )
  
  output$Abundance_thresholds <- renderUI(
    {
      if(input$which %in% c("Propagule supply")){
        textInput("abundance_thresholds", "The abundance thresholds for each proxy level", "")
      }
    }
  )
  
  output$Distance_thresholds <- renderUI(
    {
      if(input$which %in% c("Zoochory", "Hydrochory", "Agochory", "Anemochory")){
        textInput("distance_thresholds", "The distance thresholds for each proxy level", "")
      }
    }
  )
  
  output$Proxy_levels <- renderUI(
    {
      if(input$which %in% c("Propagule supply", "Zoochory", "Hydrochory", "Anemochory", "Agochory")){
        textInput("proxy_levels", "The proxy levels corresponding to the distance or abundance thresholds", "")
      }
    }
  )
  
  output$Max_radius <- renderUI(
    {
      if(input$which %in% c("Propagule supply", "Zoochory", "Hydrochory", "Anemochory", "Agochory", "Crop to extent")){
        numericInput("max_radius", "Maximum radius (m)", 1000, min = 0, max = 1e6)
      }
    }
  )
  
  output$Output_name <- renderUI(
    if(input$which != "None"){
      textInput("output_name", "Enter name of output file (no extension)", "Output_File")
    }
  )
  
  output$submit_button <- renderUI(
    {
      if(input$which != "None"){
        actionButton("submit", "Click here to begin geoprocessing")
      }
    }
  )
  
  output$download_button <- renderUI(
    {
      if(input$which != "None"){
        downloadButton("Download", "Click this button to download outputs")
      }
    }
  )
  
  the_data <- eventReactive(
    input$submit,
    {
      if(input$which == "Crop to extent"){
        
        # Ingest detection records
        files <- input$detections$datapath
        files_renamed <- gsub("/[0-9]{1}\\.", "/REPLACE\\.", files)
        for(i in 1:length(files)){
          file.rename(files[i], files_renamed[i])
        }
        detections <- shapefile(files_renamed[grep(".shp$", files_renamed)])
        
        # Ingest generic raster
        generic_raster <- raster(input$generic_raster$datapath)
        
        # Max buffer distance
        max_radius <- input$max_radius
        
        # Buffer and crop
        buffered <- gBuffer(detections, FALSE, width = max_radius)
        result <- crop(generic_raster, extent(buffered))
        
      }
      if(input$which == "Propagule supply"){
        
        # Ingest detection records
        files <- input$detections$datapath
        files_renamed <- gsub("/[0-9]{1}\\.", "/REPLACE\\.", files)
        for(i in 1:length(files)){
          file.rename(files[i], files_renamed[i])
        }
        detections <- shapefile(files_renamed[grep(".shp$", files_renamed)])
        
        # Ingest reference raster
        reference_raster <- input$reference_raster
        reference_raster <- reference_raster$datapath
        reference_raster <- raster(reference_raster)
        
        # Ingest column
        column <- input$column
        
        # Ingest distance thresholds
        abundance_thresholds <- input$abundance_thresholds
        abundance_thresholds <- as.numeric(strsplit(abundance_thresholds, ",")[[1]])
        
        # Ingest proxy levels
        proxy_levels <- input$proxy_levels
        proxy_levels <- as.numeric(strsplit(proxy_levels, ",")[[1]])
        
        # Check these
        if(any(is.na(abundance_thresholds)) | any(is.na(proxy_levels))){
          stop("Please re-enter the distance thresholds and proxy levels. Please make sure that these are entered as numbers separated by commas.")
        }
        
        # Ingest number for maximum buffer distance
        max_radius <- input$max_radius
        
        # Compute hydrochory raster
        result <- propagule_supply(detections, reference_raster, column, abundance_thresholds, proxy_levels, max_radius)
        
      }
      
      if(input$which %in% c("Zoochory", "Anemochory")){
        
        # Ingest detection records
        files <- input$detections$datapath
        files_renamed <- gsub("/[0-9]{1}\\.", "/REPLACE\\.", files)
        for(i in 1:length(files)){
          file.rename(files[i], files_renamed[i])
        }
        detections <- shapefile(files_renamed[grep(".shp$", files_renamed)])
        
        # Ingest reference raster
        reference_raster <- input$reference_raster
        reference_raster <- reference_raster$datapath
        reference_raster <- raster(reference_raster)
        
        # Ingest distance thresholds
        distance_thresholds <- input$distance_thresholds
        distance_thresholds <- as.numeric(strsplit(distance_thresholds, ",")[[1]])
        
        # Ingest proxy levels
        proxy_levels <- input$proxy_levels
        proxy_levels <- as.numeric(strsplit(proxy_levels, ",")[[1]])
        
        # Check these
        if(any(is.na(distance_thresholds)) | any(is.na(proxy_levels))){
          stop("Please re-enter the distance thresholds and proxy levels. Please make sure that these are entered as numbers separated by commas.")
        }
        
        # Ingest number for maximum buffer distance
        max_radius <- input$max_radius
        
        # Compute zoochory raster
        result <- zoochory(detections, reference_raster, distance_thresholds, proxy_levels, max_radius)
      }
      
      if(input$which %in% c("Hydrochory", "Agochory")){
        
        # Ingest stream raster
        stream_raster <- input$stream_raster
        stream_raster <- stream_raster$datapath
        stream_raster <- raster(stream_raster)
        
        # Ingest detection records
        files <- input$detections$datapath
        files_renamed <- gsub("/[0-9]{1}\\.", "/REPLACE\\.", files)
        for(i in 1:length(files)){
          file.rename(files[i], files_renamed[i])
        }
        detections <- shapefile(files_renamed[grep(".shp$", files_renamed)])
        #  distance_thresholds, proxy_levels
        
        # Ingest reference raster
        reference_raster <- input$reference_raster
        reference_raster <- reference_raster$datapath
        reference_raster <- raster(reference_raster)
        
        # Ingest distance thresholds
        distance_thresholds <- input$distance_thresholds
        distance_thresholds <- as.numeric(strsplit(distance_thresholds, ",")[[1]])
        
        # Ingest proxy levels
        proxy_levels <- input$proxy_levels
        proxy_levels <- as.numeric(strsplit(proxy_levels, ",")[[1]])
        
        # Check these
        if(any(is.na(distance_thresholds)) | any(is.na(proxy_levels))){
          stop("Please re-enter the distance thresholds and proxy levels. Please make sure that these are entered as numbers separated by commas.")
        }
        
        # Ingest number for maximum buffer distance
        max_radius <- input$max_radius
        
        # Compute hydrochory raster
        result <- hydrochory(stream_raster, detections, reference_raster, distance_thresholds, proxy_levels, max_radius)
      
      }
      
      if(input$which == "Stream edge to raster"){
        
        # Ingest shapefile
        files <- input$stream_shapefile$datapath
        files_renamed <- gsub("/[0-9]{1}\\.", "/REPLACE\\.", files)
        for(i in 1:length(files)){
          file.rename(files[i], files_renamed[i])
        }
        stream_shapefile <- shapefile(files_renamed[grep(".shp$", files_renamed)])
        
        # Ingest reference raster
        reference_raster <- raster(input$reference_raster$datapath)
        
        # Mask over reference raster
        reference_raster[] <- 1
        result <- mask(reference_raster, stream_shapefile)
        
      }
      
      if(input$which == "Recode stream raster"){
        
        # Ingest inputs
        stream_raster <- input$stream_raster
        stream_raster <- stream_raster$datapath
        stream_raster <- raster(stream_raster)
        stream_value <- input$stream_value
        # Return recoded stream raster
        result <- binary_stream_to_unary_stream(stream_raster, stream_value)
        
      }
      
      if(input$which == "Project raster"){
        
        # Load data
        project_raster <- raster(input$generic_raster$datapath)
        
        # Project the raster
        result <- projectRaster(project_raster, crs = CRS(input$proj4string))
        
      }
      
      if(input$which == "Project shapefile"){

        # Ingest shapefile
        files <- input$generic_shapefile$datapath
        files_renamed <- gsub("/[0-9]{1}\\.", "/REPLACE\\.", files)
        for(i in 1:length(files)){
          file.rename(files[i], files_renamed[i])
        }
        project_shapefile <- shapefile(files_renamed[grep(".shp$", files_renamed)])
        
        # Ingest raster with target coordinate reference system
        project_raster <- raster(input$generic_raster$datapath)

        # Project it
        result <- spTransform(project_shapefile, crs(project_raster))

      }
      
      result

    }
  )
  
  output$plot <- renderPlot(
    {
      plot(the_data())
    }
  )
  
  output$Download <- downloadHandler(
    
    filename = function(){
      if(input$which != "Project shapefile"){
        paste0(input$output_name, ".tif")
      } else {
        paste0(input$output_name, ".zip")
      }
      
    },
    content = function(file){
      if(length(Sys.glob("*.tif")) > 0){
        file.remove(Sys.glob("*.tif"))
      }
      if(input$which != "Project shapefile"){
        writeRaster(the_data(), file)
      } else {
        shapefile(the_data(), paste0(input$output_name, ".shp"), overwrite = TRUE)
        zip(zipfile = file, files = Sys.glob(paste0("*", input$output_name, ".*")))
      }
      
    }
    
  )
  
}

shinyApp(ui, server)