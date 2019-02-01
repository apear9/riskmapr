## Back end for suitability/susceptibility BN

library(shiny)
library(dplyr)
library(tidyr)
library(raster)
library(zip)
library(sp)
library(truncnorm)
library(stringr)
# library(igraph)
library(visNetwork)

## Set up option for maximum file input size
options(shiny.maxRequestSize=50*1024^2) # please, please, please try to be conservative with upload sizes

ui <- fluidPage(
  
  # App title
  
  titlePanel("Bayesian network for invasion suitability and susceptibility"),
  
  # Sidebar panel for inputs ----
  
  sidebarLayout(
    
    sidebarPanel(
      
      # Input: Select a file ----
      
      fileInput(
        "establishment", 
        "Spatial proxies for establishment:",
        multiple = TRUE,
        accept = c(".tif")
      ),
      
      helpText("Click on 'Browse...'. A new window should pop up, where you can select some raster files (with a .tif extension) to upload. Select these files and click 'Open'. Note that, if multiple files are selected, the files will be uploaded in alphabetical order by filename."),
      
      textInput(
        "establishment_weights",
        "Choose weights for establishment proxies:",
        "Enter values separated by commas or other common separators"
      ),
      
      helpText("Delete the text currently in the text input field, and replace it with a set of numbers from 1 to 3, separated by commas. These should correspond to the alphabetically ordered rasters you uploaded in the field above."),
      
      numericInput(
        "est_sd", 
        "Standard deviation for the establishment node:",
        value = 15,
        min = 0.1,
        max = 1000
      ),
      
      helpText("You can change the number above to any value between, and including, 0.1 and 1000. Don't change anything if you want to accept the default value of 15."),
      
      hr(),
      
      fileInput(
        "persistence", 
        "Spatial proxies for persistence:",
        multiple = TRUE,
        accept = c(".tif")
      ),
      
      textInput(
        "persistence_weights",
        "Choose weights for persistence proxies:",
        "Enter values separated by commas or other common separators"
      ),
      
      numericInput(
        "per_sd", 
        "Standard deviation for the persistence node:",
        value = 15,
        min = 0.1,
        max = 1000
      ),
      
      hr(),
      
      fileInput(
        "propagule_pressure", 
        "Spatial proxies for propagule pressure:",
        multiple = TRUE,
        accept = c(".tif")
      ),
      
      textInput(
        "propagule_weights",
        "Choose weights for propagule pressure proxies:",
        "Enter values separated by commas or other common separators"
      ),
      
      numericInput(
        "prg_sd", 
        "Standard deviation for the propagule pressure node:",
        value = 15,
        min = 0.1,
        max = 1000
      ),
      
      hr(),
      
      numericInput(
        "suitability_sd", 
        "Standard deviation for the suitability node:",
        value = 10,
        min = 0.1,
        max = 1000
      ),
      
      numericInput(
        "susceptibility_sd", 
        "Standard deviation for the susceptiblity node:",
        value = 10,
        min = 0.1,
        max = 1000
      ),
      
      hr(),
      
      textInput(
        "suit_name",
        "Optional: enter file name for suitability raster, no file extension",
        "Suitability"
      ),
      
      textInput(
        "susc_name",
        "Optional: enter file name for susceptibility raster, no file extension",
        "Susceptibility"
      ),
      
      helpText("You can choose the names of the downloadable files produced by this app by replacing the values in the two fields above with your own text."),
      
      hr(),
      
      actionButton("validate", "Click to view network model (validate)"),
      
      helpText("Clicking on the button above will cause the app to display a graph of the Bayesian network you have constructed by uploading files and entering inputs. Use this to check that the structure of the network is fine."),
      
      actionButton("submit", "Click to produce risk maps (run program)"),
      
      helpText("Clicking on the button above will cause the app to compute the invasion risk maps for your study area. Once this is finished, you will see a plot of the study area appear to the right."),
      
      downloadButton(outputId = "downloadData", label = "DOWNLOAD OUTPUTS"),
      
      helpText("Click on the above button to download a zip folder containing the spatial risk maps and uncertainty maps. This will only work if you can see the plot of the study area to the right.")
      
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
    req(input$propagule_pressure)
    req(input$propagule_weights)
    
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
    propagule <- input$propagule_pressure
    propagule <- propagule$name
    propagule <- gsub(".tif", "", propagule)
    propagule_wts <- input$propagule_weights
    propagule_wts <- str_split(propagule_wts, "[,/;\t]{1}")[[1]]
    propagule_wts <- as.numeric(propagule_wts)
    
    ### Lay out basic network structure (five essential nodes)
    basic_network <- data.frame(
      from = c(1, 2, 3, 4),
      to = c(3, 3, 5, 5),
      arrows = "to"
    )
    
    vertex_info <- data.frame(
      id = 1:5, 
      label = c(
        "Establishment", 
        "Persistence",
        "Suitability", 
        "Propagule\npressure", 
        "Susceptibility"
      ),
      value = 1:5,
      group = "Not user-specified",
      shape = rep("box", 5),
      color = "black",
      font.color = "white",
      shadow = rep(TRUE, 5)
    )
    
    ### Add to basic network structure based on inputs
    n_est <- length(establishment)
    n_per <- length(persistence)
    n_prg <- length(propagule)
    est_id <- 6:(5 + n_est)
    per_id <- (max(est_id) + 1):(max(est_id) + n_per)
    prg_id <- (max(per_id) + 1):(max(per_id) + n_prg)
    
    new_connections <- data.frame(
      from = c(est_id, per_id, prg_id),
      to   = c(rep(1, n_est), rep(2, n_per), rep(4, n_prg)),
      arrows = "to"
    )
    
    n_elem <- length(c(establishment, persistence, propagule))
    new_vertex_info <- data.frame(
      id = c(est_id, per_id, prg_id), 
      label = c(establishment, persistence, propagule),
      value = c(est_id, per_id, prg_id),
      group = paste("Weight =", c(establishment_wts, persistence_wts, propagule_wts)),
      shape = rep("ellipse", n_elem),
      color = colour_labeller_vectorised(c(establishment_wts, persistence_wts, propagule_wts)),
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
    req(input$propagule_pressure)
    req(input$propagule_weights)
    
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
    propagule <- input$propagule_pressure
    propagule <- propagule$datapath
    propagule_wts <- input$propagule_weights
    propagule_wts <- str_split(propagule_wts, "[,/;\t]{1}")[[1]]
    propagule_wts <- as.numeric(propagule_wts)
    propagule_sd <- input$prg_sd
    
    ## Standard deviations of suitability and susceptibility nodes
    suitability_sd <- input$suitability_sd
    susceptibility_sd <- input$susceptibility_sd
    
    ## Find length of names
    nn_per <- length(persistence)
    nn_est <- length(establishment)
    nn_pgl <- length(propagule)
    
    ## Check that lengths are what they should be
    
    if(nn_per != length(persistence_wts)){
      stop("The number of persistence weights is not equal to the number of proxy rasters provided.")
    }
    if(nn_est != length(establishment_wts)){
      stop("The number of establishment weights is not equal to the number of proxy rasters provided.")
    }
    if(nn_pgl != length(propagule_wts)){
      stop("The number of propagule weights is not equal to the number of proxy rasters provided.")
    }
    
    ## Construct indices, persistence first
    i_per <- 1:nn_per
    i_est <- (nn_per + 1):(nn_per + nn_est)
    i_prg <- (nn_per + nn_est + 1):(nn_per + nn_est + nn_pgl)
    
    ## Read in rasters as stack
    drop_ind <- 1:(length(c(persistence, establishment, propagule)) - 1)
    suit_ras <- stack(c(persistence, establishment, propagule))
    suit_ras_df <- as.data.frame(suit_ras)
    suit_ras <- dropLayer(suit_ras, drop_ind)
    
    ## Extract distinct rows WITHOUT rows involving NAs.
    # suit_ras_df_dn <- suit_ras_df[!duplicated(suit_ras_df), ]
    suit_ras_df_dn <- distinct(suit_ras_df)
    # suit_ras_df_dn <- suit_ras_df_dn[apply(suit_ras_df_dn, 1, function(x) !any(is.na(x))), ]
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
    # prg_cols <- as.matrix(suit_ras_df_dn[, i_prg])
    prg_wets <- propagule_wts
    
    ## Compute distribution of suitability 
    
    # Empty numeric vectors
    st <- st_sd <- sc <- sc_sd <- numeric(nrow(suit_ras_df_dn))
    
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
      
      # Now construct the propagule pressure node
      # prg_vars <- prg_cols[i, ]
      prg_vars <- suit_ras_df_dn[i, i_prg]
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
    
    # Derive ID column
    # suit_ras_df_dn$id <- apply(suit_ras_df_dn, 1, function(x) paste(x, collapse = ""))
    suit_ras_df_dn <- unite(suit_ras_df_dn, "id", sep = "", remove = TRUE)
    suit_ras_df_dn$Suitability <- st
    suit_ras_df_dn$Suitability_SD <- st_sd
    suit_ras_df_dn$Susceptibility <- sc
    suit_ras_df_dn$Susceptibility_SD <- sc_sd
    
    # Join back to full dataset
    # suit_ras_df$id <- apply(suit_ras_df, 1, function(x) paste(x, collapse = ""))
    suit_ras_df <- unite(suit_ras_df, "id", sep = "", remove = TRUE)
    s_result <- left_join(suit_ras_df, suit_ras_df_dn, by = "id")
    
    rm(suit_ras_df, suit_ras_df_dn)
    
    ### Put back into raster
    suit_ras$Suitability <- s_result$Suitability
    suit_ras$Suitability_SD <- s_result$Suitability_SD
    suit_ras$Susceptibility <- s_result$Susceptibility
    suit_ras$Susceptibility_SD <- s_result$Susceptibility_SD
    
    rm(s_result, st, st_sd, sc, sc_sd)
    
    suit_ras
    
  })
  
  output$mainplot <- renderPlot(
    {
      
      ### Plot
      spplot(the_plots(), c("Suitability", "Susceptibility"))
      
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
      susc_fn <- paste0(input$susc_name, ".tif")
      susc_sd_fn <- paste0(input$susc_name, "_SD.tif")
      writeRaster(the_stack$Suitability, suit_fn)
      writeRaster(the_stack$Susceptibility, susc_fn)
      writeRaster(the_stack$Suitability_SD, suit_sd_fn)
      writeRaster(the_stack$Susceptibility_SD, susc_sd_fn)
      zip(zipfile = file, files = Sys.glob("*.tif"))
    },
    contentType = "application/zip"
    
  )
  
}

shinyApp(ui, server)