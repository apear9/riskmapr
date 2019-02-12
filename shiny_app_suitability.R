## Back end for suitability/susceptibility BN

library(rgdal)
library(shiny)
library(raster)
library(zip)
library(sp)
library(truncnorm)
library(stringr)
library(visNetwork)

## Set up option for maximum file input size
options(shiny.maxRequestSize=50*1024^2) # please, please, please try to be conservative with upload sizes

ui <- fluidPage(
  
  # App title
  
  titlePanel("Rapid weed riskmapr (suitability model)"),
  
  # Sidebar panel for inputs ----
  
  sidebarLayout(
    
    sidebarPanel(width = 6,
                 
                 # Input: Select a file ----
                 
                 fileInput(
                   "establishment", 
                   "Spatial proxies for risk factors (establishment)",
                   multiple = TRUE,
                   accept = c(".tif")
                 ),
                 
                 helpText("Upload pre-processed spatial proxies for all identified risk factors affecting plant establishment. Select all relevant files (must have .TIF extension) at once and click 'Open'. Files are automatically uploaded in alphabetical order. Upload limit is 50MB, but functionality has only been confirmed for total upload sizes <20MB."),
                 
                 textInput(
                   "establishment_weights",
                   "Risk factor weights (establishment)"
                 ),
                 
                 helpText("Enter weights for all identified risk factors affecting plant establishment. Weights must equal '1', '2' or '3', be separated by commas and ordered alphabetically by spatial proxy name."),
                 
                 numericInput(
                   "est_sd", 
                   "Standard deviation (establishment)",
                   value = 15,
                   min = 0.1,
                   max = 100
                 ),
                 
                 helpText("Enter the standard deviation used for computing the probability distribution of plant establishment as a function of its weighted risk factors. The default is '15'. This can be changed to any reasonable value, keeping in mind that the mean is between 0 and 100 (depending on the state of each risk factor)"),
                 
                 fileInput(
                   "persistence", 
                   "Spatial proxies for risk factors (persistence)",
                   multiple = TRUE,
                   accept = c(".tif")
                 ),
                 
                 helpText("Upload pre-processed spatial proxies for all identified risk factors affecting plant persistence. For details, see above."),
                 
                 textInput(
                   "persistence_weights",
                   "Risk factor weights (persistence)"
                 ),
                 
                 helpText("Enter the weights for all identified risk factors affecting plant persistence. For details, see above."),
                 
                 numericInput(
                   "per_sd", 
                   "Standard deviation (persistence)",
                   value = 15,
                   min = 0.1,
                   max = 100
                 ),
                 
                 helpText("Enter the standard deviation used for computing the probability distribution of plant persistence. For details, see above."),
                 
                 numericInput(
                   "suitability_sd", 
                   "Standard deviation (suitability)",
                   value = 10,
                   min = 0.1,
                   max = 1000
                 ),
                 
                 helpText("Enter the standard deviation used for computing the probability distribution of invasion risk (suitability) as a function of plant establishment and persistence. The default is '10'. This is lower than SD = '15' above in order to limit the propagated uncertainty in the model, but can be changed to any reasonable value."),
                 
                 textInput(
                   "suit_name",
                   "Optional: name risk map (suitability)",
                   "Suitability"
                 ),
                 
                 helpText("Choose a descriptive name for the generated risk map (no file extension). Please choose the file name before running the tool."),
                 
                 actionButton("validate", "VISUALIZE RISK MODEL"),
                 
                 helpText("Click to visualize and validate the structure of your risk model (suitability). The model is displayed on the right-hand panel, showing uploaded spatial proxies colour-coded by assigned risk factor weights."),
                 
                 actionButton("submit", "RUN RISK MODEL"),
                 
                 helpText("Click to run your risk model (suitability). Two spatial files (.TIF) are generated: a suitability index map (the model expected value), and an uncertainty map (the model standard deviation) This may take several minutes, depending on the size of spatial proxies. Once completed, the risk map is displayed on the right-hand panel."),
                 
                 downloadButton(outputId = "downloadData", label = "DOWNLOAD RISK MAP"),
                 
                 helpText("Once the risk map has been generated and displayed, click to download a .ZIP folder with model outputs (suitability index map + uncertainty map).")
                 
    ),
    
    mainPanel(width = 6, 
              visNetworkOutput("valiplot"),
              plotOutput("mainplot")
    )
    
  )
  
)


server <- function(input, output){
  
  the_graph <- eventReactive(input$validate, {
    
    ### Two functions needed for colourmapping the network edges
    colour_labeller <- function(wt) switch(
      wt, "1" = "forestgreen", "2" = "orange", "3" = "red"
    )
    colour_labeller_vectorised <- function(wts){
      if(!all(wts %in% 1:3)){
        stop("All weights must be 1, 2, or 3")
      }
      unlist(
        lapply(wts, colour_labeller)
      )
    }
    
    ### Get inputs for model validation
    req(input$persistence)
    req(input$persistence_weights)
    req(input$establishment)
    req(input$establishment_weights)
    
    ### Preprocess the inputs 
    persistence <- input$persistence
    persistence <- persistence$name
    persistence <- gsub(".tif", "", persistence)
    persistence_wts <- input$persistence_weights
    persistence_wts <- str_split(persistence_wts, "[,/;\t]{1}")[[1]]
    persistence_wts <- as.numeric(persistence_wts)
    establishment <- input$establishment
    establishment <- establishment$name
    establishment <- gsub(".tif", "", establishment)
    establishment_wts <- input$establishment_weights
    establishment_wts <- str_split(establishment_wts, "[,/;\t]{1}")[[1]]
    establishment_wts <- as.numeric(establishment_wts)
    
    ### Lay out basic network structure (five essential nodes)
    basic_network <- data.frame(
      from = c(1, 2),
      to = c(3, 3),
      arrows = "to"
    )
    vertex_info <- data.frame(
      id = 1:3, 
      label = c(
        "Establishment", 
        "Persistence",
        "Suitability"
      ),
      value = 1:3,
      group = "Not user-specified",
      shape = rep("box", 3),
      color = "black",
      font.color = "white",
      shadow = rep(TRUE, 3)
    )
    
    ### Add to basic network structure based on inputs
    n_est <- length(establishment)
    n_per <- length(persistence)
    est_id <- 4:(3 + n_est)
    per_id <- (max(est_id) + 1):(max(est_id) + n_per)
    
    new_connections <- data.frame(
      from = c(est_id, per_id),
      to   = c(rep(1, n_est), rep(2, n_per)),
      arrows = "to"
    )
    
    n_elem <- length(c(establishment, persistence))
    new_vertex_info <- data.frame(
      id = c(est_id, per_id), 
      label = c(establishment, persistence),
      value = c(est_id, per_id),
      group = paste("Weight =", c(establishment_wts, persistence_wts)),
      shape = rep("ellipse", n_elem),
      color = colour_labeller_vectorised(c(establishment_wts, persistence_wts)),
      font.color = "white",
      shadow = rep(FALSE, n_elem)
    )
    
    all_network <- rbind(
      basic_network,
      new_connections
    )
    all_vertex <- rbind(
      vertex_info,
      new_vertex_info
    )
    the_graph <- visNetwork(all_vertex, all_network, height = "600px", width = "100%") %>% #,
      visGroups(groupname = "Weight = 1", color = "forestgreen", font = list(color = "white")) %>%
      visGroups(groupname = "Weight = 2", color = "orange", font = list(color = "white")) %>%
      visGroups(groupname = "Weight = 3", color = "red", font = list(color = "white")) %>%
      visGroups(groupname = "Not user-specified", color = "black", shape = "box", font = list(color = "white")) %>%
      visLegend(main = "Legend") %>%
      visPhysics(enabled = FALSE)
    
    the_graph
    
  }
  )
  
  output$valiplot <- renderVisNetwork(
    {
      
      the_graph()
      
    }
  )
  
  the_plots <- eventReactive(input$submit, {
    
    ### Define functions for finding expectations and standard deviations
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
    } # Standard definitions for discrete random variables
    
    # Stripped down unique function from raster that does not allow for in-memory processing
    unique_out_of_memory <- function(x){
      # MODIFIED SOURCE CODE FROM THE PACKAGE 'RASTER', FROM FUNCTION raster::unique().
      # MODIFIED 5 FEB, 2019
      nl <- nlayers(x)
      un <- list(length = nl, mode = "list")
      tr <- blockSize(x, n = nl, minblocks = nl * 10)
      un <- NULL
      for (i in 1:tr$n) {
        v <- dplyr::distinct(
          as.data.frame(
            getValues(x, row = tr$row[i], nrows = tr$nrows[i])
          )
        )
        un <- rbind(v, un)
      }
      return(un)
    }
    
    ### Get inputs
    req(input$persistence)
    req(input$persistence_weights)
    req(input$establishment)
    req(input$establishment_weights)
    
    ### Some preprocessing to force the files to be loaded in alphanumeric order
    per_an <- input$persistence
    per_an <- per_an$name
    per_an <- gsub(".tif", "", per_an) # Not that it really matters...
    est_an <- input$establishment
    est_an <- est_an$name
    est_an <- gsub(".tif", "", est_an)
    per_or <- order(per_an)
    est_or <- order(est_an)
    
    ### Get persistence and establishment as rasters
    persistence <- input$persistence
    persistence <- persistence$datapath[per_or]
    persistence_wts <- input$persistence_weights
    persistence_wts <- str_split(persistence_wts, "[,/;\t]{1}")[[1]]
    persistence_wts <- as.numeric(persistence_wts)
    persistence_sd <- input$per_sd
    establishment <- input$establishment
    establishment <- establishment$datapath[est_or]
    establishment_wts <- input$establishment_weights
    establishment_wts <- str_split(establishment_wts, "[,/;\t]{1}")[[1]]
    establishment_wts <- as.numeric(establishment_wts)
    establishment_sd <- input$est_sd
    
    ## Standard deviations of suitability and susceptibility nodes
    suitability_sd <- input$suitability_sd
    
    ## Find length of names
    nn_per <- length(persistence)
    nn_est <- length(establishment)
    
    ## Check that lengths are what they should be
    if(nn_per != length(persistence_wts)){
      stop("The number of persistence weights is not equal to the number of proxy rasters provided.")
    }
    if(nn_est != length(establishment_wts)){
      stop("The number of establishment weights is not equal to the number of proxy rasters provided.")
    }
    
    ## Construct indices, persistence first
    i_per <- 1:nn_per
    i_est <- (nn_per + 1):(nn_per + nn_est)
    
    ## Read in rasters as stack
    suit_ras <- stack(c(persistence, establishment))
    
    # Find unique combinations of values
    message("Finding the unique combinations of proxies")
    suit_ras_df_dn <- unique_out_of_memory(suit_ras) # ... can take a while, but doesn't break the bank when it comes to memory usage.
    suit_ras_df_dn <- dplyr::distinct(suit_ras_df_dn)
    
    # Remove any rows with NAs or values outside the 0, 100 range
    message("Getting rid of meaningless combinations of values")
    ind_na <- rowSums(is.na(suit_ras_df_dn)) == 0
    suit_ras_df_dn <- suit_ras_df_dn[ind_na, ]
    ind_rn <- rowSums(suit_ras_df_dn < 0 | suit_ras_df_dn > 100) == 0
    suit_ras_df_dn <- suit_ras_df_dn[ind_rn, ]
    rm(ind_rn, ind_na)
    gc()
    
    # Subsets as required for the analysis
    per_wets <- persistence_wts 
    est_wets <- establishment_wts 
    
    # Empty numeric vectors, needed for loop
    st <- st_sd <- numeric(nrow(suit_ras_df_dn))
    
    # Main loop
    message("Starting the main loop")
    for(i in 1:nrow(suit_ras_df_dn)){
      
      # Establishment
      est_vars <- suit_ras_df_dn[i, i_est]
      est_mean <- sum(est_vars * est_wets)/sum(est_wets)
      est <- dtruncnorm(seq(0, 100, 25), 0, 100, est_mean, establishment_sd)
      est <- est/sum(est)
      names(est) <- seq(0, 100, 25)
      
      # Persistence
      per_vars <- suit_ras_df_dn[i, i_per]
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
      
      # Take expectation as prediction
      st[i] <- exp_discrete(suit)
      st_sd[i] <- std_discrete(suit)
      
    }
    
    # Derive ID column
    message("Exited from main loop")
    suit_ras_df_dn <- as.data.frame(suit_ras_df_dn)
    suit_ras_df_dn$Suitability <- st 
    rm(st)
    suit_ras_df_dn$Suitability_SD <- st_sd
    rm(st_sd)
    gc()
    
    # Begin the process of joining this back to the full dataset, all done by manipulating the files and without ingesting the entire raster into memory
    chunk_info <- blockSize(suit_ras, n = nlayers(suit_ras), minblocks = nlayers(suit_ras) * 10) # the n_layers * 10 thing is arbitrary, just trying to make sure R doesn't bite off more than it can chew.
    
    # Prepare to write by constructing file names
    suit_fn <- paste0(input$suit_name, ".tif")
    suit_sd_fn <- paste0(input$suit_name, "_SD.tif")
    
    # Remove existing .tif files so they don't get packaged up 
    if(length(Sys.glob("*.tif")) > 0){
      file.remove(Sys.glob("*.tif"))
    }
    
    # Open file connections
    message("Preparing to write rasters for suitability and uncertainty.")
    f1 <- writeStart(suit_ras[[1]], suit_fn, overwrite = TRUE)
    f2 <- writeStart(suit_ras[[1]], suit_sd_fn, overwrite = TRUE)
    
    # Then the loop, ingesting the raster by chunks, writing it by the same chunks
    for(i in 1:chunk_info$n){
      tmp_df <- as.data.frame(
        getValues(suit_ras, row = chunk_info$row[i], nrows = chunk_info$nrows[i])
      )
      vals_df <- dplyr::left_join(
        tmp_df, 
        suit_ras_df_dn, 
        by = names(tmp_df)
      )
      rm(tmp_df)
      gc()
      f1 <- writeValues(f1, vals_df$Suitability, chunk_info$row[i])
      f2 <- writeValues(f2, vals_df$Suitability_SD, chunk_info$row[i])
      rm(vals_df)
      gc()
      
    }
    f1 <- writeStop(f1)
    f2 <- writeStop(f2)
    rm(f1, f2)
    gc()
    
    ### Put back into raster
    message("Raster ready for display")
    suit_ras <- raster(suit_fn)
    suit_ras
    
  })
  
  output$mainplot <- renderPlot(
    {
      
      ### Plot
      spplot(the_plots())
      
    }
  )
  
  output$downloadData <- downloadHandler(
    
    filename = "Raster_Exports.zip",
    content = function(file){
      # Since all the files have already been created, all we have to do is zip them up.
      zip(zipfile = file, files = Sys.glob(paste0(input$suit_name, "*")))
    },
    contentType = "application/zip"
    
  )
  
}

shinyApp(ui, server)