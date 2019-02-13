# The riskmapr shiny applications

This is the GitHub repository for the suite of riskmapr shiny apps.

There are three .R files corresponding to the three Shiny apps:

* shiny_app_geoprocessing.R, which provides the common geoprocessing functions required to prepare the data for use in the other apps.
* shiny_app_suitability.R, which allows users to create weed invasion suitability maps and download the rasters. 
* shiny_app_susceptibility.R, which allows users to create weed invasion suitability and susceptibility maps and download the rasters.

There are two additional scripts which can be run to ensure the R packages required for the apps are installed:

* installing_packages_for_geoprocessing.R (for shiny_app_geoprocessing.R)
* installing_packages_for_suitability_susceptibility.R (for shiny_app_suitability.R and shiny_app_susceptibility.R)

# Running the shiny apps

The process for running these Shiny apps is easy if you are using RStudio. Simply:

* Download the code for the apps.
* Install the requisite packages using the installing_packages scripts.
* Open the .R file for the app you want to run.
* Go to the 'Run App' button in the top-right of the script window.
* Click on the small black arrow to the right of the text which says 'Run App'.
* Click 'Run External' from the small drop-down menu that appears. This will allow the app to run in an internet browser, which in turn will allow the outputs to be downloaded without fuss. 
* Click 'Run App'.

