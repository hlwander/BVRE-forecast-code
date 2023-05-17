# BVRE-forecast-code

This branch includes code required to reproduce figures from BVR FLARE DA experiments in Wander et al. *Ecosphere* paper titled "Data assimilation experiments inform monitoring needs for near-term ecological forecasts in a eutrophic reservoir"

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
