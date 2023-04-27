# BVRE-forecast-code

This branch includes code required to reproduce figures from BVR FLARE DA experiments in Wander et al. *Ecosphere* paper titled "Data assimilation experiments inform monitoring needs for near-term, ecological forecasts in a eutrophic reservoir"

# Instructions to reproduce manuscript + SI figures:

1.  Download or clone github repository to your local computer
2.  Download forecast and score parquet files generated from combined_workflow.R script in the workflows/DA_experiments folder **(LINK TO OTHER ZENODO PUB?**)
3.  Run "BVR_FLARE_ms_figs.R" script in the workflows/DA_experiments folder to reproduce manuscript figures
4.  Run "BVR_FLARE_UC_figs.R" script in the workflows/DA_experiments folder to reproduce Fig. 9 (proportion of IC uncertainty) and SI figures for forecasts run without initial conditions uncertainty
5.  Run "BVR_FLARE_SI_figs.R" script in the workflows/DA_experiments folder to reproduce remaining SI figures

# Instructions to reproduce FLARE forecasts:

1.   Run "install.R" in the workflows/DA_experiments folder to download GLM and FLARE packages and their dependencies

2.   Run "01_generate_targets.R" in the workflows/DA_experiments folder to download driver data and observation files from EDI or GitHub

3.  Run "02_combined_workflows.R" in the workflows/DA_experiments folder to iteratively generate forecasts for every data assimilation frequency and day in the forecast period

    **Note** Running forecasts for 365 days and all four data assimilation frequencies will take \> 10 days and the forecasts that are generated will be slightly different than the ones used in this analysis due to the stochasticity involved in accounting for uncertainty in the FLARE workflow.
