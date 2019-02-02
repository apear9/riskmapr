## Back end for suitability/susceptibility BN

library(shiny)
library(dplyr)
library(tidyr)
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
  
  titlePanel("Rapid riskmapr - suitability for weed invasion"),
  
  # Sidebar panel for inputs ----
  
  sidebarLayout(
    
    sidebarPanel(
      
      # Input: Select a file ----
      
      fileInput(
        "establishment", 
        "Spatial proxies for risk factors (establishment)",
        multiple = TRUE,
        accept = c(".tif")
      ),

      helpText("Upload pre-processed spatial proxies for all identified risk factors affecting plant establishment. Select all relevant files (must have .TIF extension) at once and click 'Open'. Files are automatically uploaded in alphabetical order. Upload limit is 50MB, but functionality has only been confirmed for <20MB."),
      
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
        max = 1000
      ),
      
      helpText("Enter the standard deviation used for computing the probability distribution of plant establishment as a function of its weighted risk factors. The default is '15'. This can be changed to any reasonable value, keeping in mind that the mean is between 0 and 100 (depending on the state of each risk factor)"),
      
      hr(),
      
      fileInput(
        "persistence", 
        "Spatial proxies for risk factors (persistence)",
        multiple = TRUE,
        accept = c(".tif")
      ),
      
      helpText("Upload pre-processed spatial proxies for all identified risk factors affecting plant persistence. Details see above."),
      
      textInput(
        "persistence_weights",
        "Risk factor weights (persistence)"
      ),
      
      helpText("Enter the weights for all identified risk factors affecting plant persistence. Details see above"),
      
      numericInput(
        "per_sd", 
        "Standard deviation (persistence)",
        value = 15,
        min = 0.1,
        max = 1000
      ),
      
      helpText("Enter the standard deviation used for computing the probability distribution of plant persistence. Details see above."),
      
      hr(),
      
      numericInput(
        "suitability_sd", 
        "Standard deviation (suitability)",
        value = 10,
        min = 0.1,
        max = 1000
      ),
      
      helpText("Enter the standard deviation used for computing the probability distribution of invasion risk (suitability) as a function of plant establishment and persistence. The default is '10'. This is lower than SD = '15' above in order to limit the propagated uncertainty in the model, but can be changed to any reasonable value."),
      
      hr(),
      
      actionButton("validate", "Visualize risk model"),
      
      helpText("Click to visualize and validate the structure of your risk model (suitability). The model is displayed on the right-hand panel, showing uploaded spatial proxies and risk factor weights (colour-coded network links)."),
      
      actionButton("submit", "Run risk model"),
      
      helpText("Click to run your risk model (suitability). Two spatial files (.TIF) are generated: a suitability index map (the model expected value), and an uncertainty map (the model standard deviation) This may take several minutes, depending on the size of spatial proxies. Once completed, the risk map is displayed on the right-hand panel."),
      
      hr(),
      
      textInput(
        "suit_name",
        "Optional: name risk map (suitability)",
        "Suitability"
      ),
      
      helpText("Choose a descriptive name for the generated risk map before downloading (no file extension)."),
      
      hr(),
      
      downloadButton(outputId = "downloadData", label = "Download risk map"),
      
      helpText("Once the risk map has been generated and displayed, click to download a .ZIP folder with model outputs (suitability index map + uncertainty map).")
      
    ),
    
    mainPanel(
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
    }
    
    ### Get inputs
    req(input$persistence)
    req(input$persistence_weights)
    req(input$establishment)
    req(input$establishment_weights)
    
    ### Get persistence and establishment as rasters
    persistence <- input$persistence
    persistence <- persistence$datapath
    persistence_wts <- input$persistence_weights
    persistence_wts <- str_split(persistence_wts, "[,/;\t]{1}")[[1]]
    persistence_wts <- as.numeric(persistence_wts)
    persistence_sd <- input$per_sd
    establishment <- input$establishment
    establishment <- establishment$datapath
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
    drop_ind <- 1:(length(c(persistence, establishment)) - 1)
    suit_ras <- stack(c(persistence, establishment))
    suit_ras_df <- as.data.frame(suit_ras)
    suit_ras <- dropLayer(suit_ras, drop_ind)
    
    ## Extract distinct rows WITHOUT rows involving NAs.
    #suit_ras_df_dn <- suit_ras_df[!duplicated(suit_ras_df), ]
    suit_ras_df_dn <- distinct(suit_ras_df)
    # suit_ras_df_dn <- na.omit(suit_ras_df_dn)[apply(suit_ras_df_dn, 1, function(x) !any(is.na(x))), ]
    suit_ras_df_dn <- na.omit(suit_ras_df_dn)
    # suit_ras_df_dn <- suit_ras_df_dn[apply(
    #   suit_ras_df_dn, 
    #   1, 
    #   function(x) !(any(x > 100) | any(x < 0))), ]
    suit_ras_df_dn <- filter_all(suit_ras_df_dn, all_vars(. >= 0 & . <= 100))
    
    # Subsets as required for the analysis
    # per_cols <- as.matrix(suit_ras_df_dn[, i_per])
    per_wets <- persistence_wts 
    # est_cols <- as.matrix(suit_ras_df_dn[, i_est])
    est_wets <- establishment_wts 
    
    ## Compute distribution of suitability 
    
    # Empty numeric vectors
    st <- st_sd <- numeric(nrow(suit_ras_df_dn))
    
    # Main loop
    for(i in 1:nrow(suit_ras_df_dn)){
      
      # Establishment
      # est_vars <- est_cols[i, ]
      est_vars <- suit_ras_df_dn[i, i_est]
      est_mean <- sum(est_vars * est_wets)/sum(est_wets)
      est <- dtruncnorm(seq(0, 100, 25), 0, 100, est_mean, establishment_sd)
      est <- est/sum(est)
      names(est) <- seq(0, 100, 25)
      
      # Persistence
      # per_vars <- per_cols[i, ]
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
      suit_wets <- sum(est_wets) + sum(per_wets)
      
      # Take expectation as prediction
      st[i] <- exp_discrete(suit)
      st_sd[i] <- std_discrete(suit)
      
    }
    
    # Derive ID column
    # suit_ras_df_dn$id <- apply(suit_ras_df_dn, 1, function(x) paste(x, collapse = ""))
    suit_ras_df_dn <- unite(suit_ras_df_dn, "id", sep = "", remove = TRUE)
    suit_ras_df_dn$Suitability <- st
    suit_ras_df_dn$Suitability_SD <- st_sd
    
    # Join back to full dataset
    # suit_ras_df$id <- apply(suit_ras_df, 1, function(x) paste(x, collapse = ""))
    suit_ras_df <- unite(suit_ras_df, "id", sep = "", remove = TRUE)
    s_result <- left_join(suit_ras_df, suit_ras_df_dn, by = "id")
    
    rm(suit_ras_df, suit_ras_df_dn)
    
    ### Put back into raster
    suit_ras$Suitability <- s_result$Suitability
    suit_ras$Suitability_SD <- s_result$Suitability_SD
    
    rm(s_result, st, st_sd)
    
    suit_ras
    
  })
  
  output$mainplot <- renderPlot(
    {

      ### Plot
      spplot(the_plots(), "Suitability")
      
    }
  )
  
  output$downloadData <- downloadHandler(
    
    filename = "Raster_Exports.zip",
    content = function(file){
      the_stack <- stack(the_plots())
      if(length(Sys.glob("*.tif")) > 0){
        file.remove(Sys.glob("*.tif"))
      }
      suit_fn <- paste0(input$suit_name, ".tif")
      suit_sd_fn <- paste0(input$suit_name, "_SD.tif")
      writeRaster(the_stack$Suitability, suit_fn)
      writeRaster(the_stack$Suitability_SD, suit_sd_fn)
      zip(zipfile = file, files = Sys.glob("*.tif"))
    },
    contentType = "application/zip"
    
  )
  
}

shinyApp(ui, server)