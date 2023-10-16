# BVRE-forecast-code       

This code reproduces figures from the Beaverdam Reservoir forecasting data assimilation experiments using the FLARE (Forecasting Lake And Reservoir Ecosystems) system in the manuscript by Wander et al. titled "Data assimilation experiments inform monitoring needs for near-term ecological forecasts in a eutrophic reservoir." If you have any questions, contact Heather Wander at hwander\@vt.edu. 

# Instructions to reproduce manuscript + SI figures:         

1.  Download or clone github repository to your local computer
2.  Run `install.R` in the `workflows/DA_experiments` folder to download GLM and FLARE packages and their dependencies
3.  Run `download_scores.R` in the `workflows/DA_experiments` folder to download driver data and observation files from EDI or GitHub
4.  Run `BVR_FLARE_ms_figs.R` script in the `workflows/DA_experiments` folder to reproduce manuscript and supplemental figures
5.  Run `BVR_FLARE_UC_figs.R` script in the `workflows/DA_experiments` folder to reproduce Fig. 9 (proportion of IC uncertainty) and SI figures for forecasts run without initial conditions uncertainty

# Instructions to reproduce FLARE forecasts and scores: 

1.  Run `install.R` in the `workflows/DA_experiments` folder to download GLM and FLARE packages and their dependencies

2.  Run `combined_workflows.R` in the `workflows/DA_experiments` folder to iteratively generate forecasts for every data assimilation frequency and day in the forecast period

    **Note** Running forecasts for 365 days and all four data assimilation frequencies will take \> 10 days.

# Instructions for reproducing using Docker    

1.  Download and install Docker to your computer (<https://www.docker.com>)

2.  At the command line, run `docker run --rm -ti -e PASSWORD=yourpassword -p 8787:8787 rqthomas/wander_et_al:latest`

3.  Open a webbrowser and enter `http://localhost:8787`. You will see an Rstudio login screen. The user name is `rstudio` and the password is `yourpassword`

4.  In the Rstudion session: File -\> Open project -\> select BVRE-forecast-code/BVRE-forecast-code.Rproj

5.  Follow the instructions above for reproducing the figures or the forecasts (**note: the R packages are already installed in the Docker container so `install.R` does not need to be run**)
