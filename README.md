# The riskmapr shiny apps

This website/GitHub repository is for the suite of riskmapr shiny apps. In this document, we outline how to download and run them. 

## Contents of the repository

The GitHub repository apear9/riskmapr contains three .R files corresponding to the three riskmapr Shiny apps:

* shiny_app_geoprocessing.R, which lets users perform a range of customized geoprocessing functions to prepare spatial data for use in the app shiny_app_susceptibility, and download the outputs.
* shiny_app_suitability.R, which lets users run a spatially explicit model of suitability for weed invasion, and download the generated risk maps. 
* shiny_app_susceptibility.R, which lets users run a spatially explicit model of susceptibility to weed invasion, and download the generated risk maps.

There are two additional scripts which can be run to ensure the R packages required for the apps are installed:

* installing_packages_for_geoprocessing.R (for shiny_app_geoprocessing.R)
* installing_packages_for_suitability_susceptibility.R (for shiny_app_suitability.R and shiny_app_susceptibility.R)

The other files in the repository exist for the purposes of licensing this software under GPL v3.0 and supporting the website that holds this documentation. These are not used to run the riskmapr apps. 

## Deploying and running the shiny apps locally

The process for running these Shiny apps is easy if you are using RStudio. Simply:

* Download the code for the apps. Please make sure the R scripts have been saved in a new folder.  

![download](https://user-images.githubusercontent.com/17267197/52981526-f4fd7480-342b-11e9-8ea4-d2e296418c6c.png)

* Install the requisite packages using the installing_packages scripts.
* Open the .R file for the app you want to run.
* Go to the 'Run App' button in the top-right of the script window. 
    * Click on the small black arrow next to the button.
    * Click 'Run External' from the small drop-down menu that appears. This will allow the app to run in an internet browser, which in turn will allow the outputs to be downloaded without fuss. 

![runexternalfurtherinstructions](https://user-images.githubusercontent.com/17267197/52981533-fd55af80-342b-11e9-82d7-374203bc4371.png)

* Click 'Run App'.

![shiny_app_run_instructions](https://user-images.githubusercontent.com/17267197/52686165-38716200-2f98-11e9-89e6-4e3e1e0f4b29.png)

## Deploying the shiny apps online to your own shinyapps.io account

There may be situations where it is desirable to share the app with collaborators without forcing them to interact with the code or the R language. In this case, one option is to deploy the app to [shinyapps.io](https://www.shinyapps.io/) (or some other R Shiny server). An account with [shinyapps.io](https://www.shinyapps.io) with space for 5 apps is free, so we provide instructions for deploying the app to this service. Go to the website and sign up following the instructions on the page.

Once the shinyapps.io account has been created, simply:

* STEP 1
* STEP 2
* STEP 3

Instructions to come. 
