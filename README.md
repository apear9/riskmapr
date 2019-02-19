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

* Make sure all required packages have been installed on your machine.
* Open the R script for the desired app in RStudio and click the blue publishing icon in the top right-hand corner of the scripting window.

![publishapp](https://user-images.githubusercontent.com/17267197/53022457-f029d680-34a6-11e9-9f85-10bf21d38de7.png)

* Click through the next few windows asking you to choose a publishing service. Choose shinyapps.io when asked which one to use. 

![connectaccount2](https://user-images.githubusercontent.com/17267197/53022651-5878b800-34a7-11e9-9191-b18646988f8b.PNG)

* You will be asked to link your shinyapps.io account to your computer so that the app can be published through RStudio. 
* You will see this screen. From here, go to an internet browser and log into shinyapps.io using your account details.

![linkshinyappsioaccount](https://user-images.githubusercontent.com/17267197/53022687-6af2f180-34a7-11e9-9bce-76a35e8a685b.PNG)

* Find the menu for the access tokens. The directions are indicated below.
    
![findaccesstokens](https://user-images.githubusercontent.com/17267197/53022855-ca510180-34a7-11e9-90da-bea95dc8d429.png)

* You will find a table with entries that look like this. Select your access token and click on the 'Show' button. You will then be taken to a pop-up window where you can click a button reading 'Show Secret'. Click on it.
    
![revealaccesstoken](https://user-images.githubusercontent.com/17267197/53022905-e654a300-34a7-11e9-8a29-0e19603c8dd7.png)

* Copy and paste the codechunk starting with `rsconnect::` into this window which you should have open in RStudio.
    
![linkshinyappsioaccount](https://user-images.githubusercontent.com/17267197/53023259-91655c80-34a8-11e9-9c3c-e6991d5e1f8c.PNG)

* Click 'Connect Account'.
* Now the publishing menu should be open. 

![publishingmenu2](https://user-images.githubusercontent.com/17267197/53023364-cd98bd00-34a8-11e9-822b-28062cd24178.png)

* Before proceeding any further, check that the app is being uploaded to the right account.
* Then, firstly, click the button 'Uncheck All' (indicated as 1. in the figure above). 
* Secondly, rename the app to something of your choosing (in the box indicated as 2. in the figure above).
* Thirdly, click the Publish button (3. in the figure above).
* R may prompt you to install packages. Do this if required. 
* A new window should open detailing the progress of RStudio in deploying the app to your shinyapps.io account. An internet browser window will open once this is complete, and you will be able to use the app online. The URL can be shared with anyone. 
