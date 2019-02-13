# The riskmapr shiny applications

This is the GitHub repository for the suite of riskmapr shiny apps.

There are three .R files corresponding to the three Shiny apps:

* shiny_app_geoprocessing.R, which lets users perform a range of customized geoprocessing functions to prepare spatial data for use in the app shiny_app_susceptibility, and download the outputs.
* shiny_app_suitability.R, which lets users run a spatially explicit model of suitability for weed invasion, and download the generated risk maps. 
* shiny_app_susceptibility.R, which lets users run a spatially explicit model of susceptibility to weed invasion, and download the generated risk maps.

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

![shiny_app_run_instructions](https://user-images.githubusercontent.com/17267197/52686165-38716200-2f98-11e9-89e6-4e3e1e0f4b29.png)
