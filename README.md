# BVRE-forecast-code

This branch includes code required to reproduce figures from BVR FLARE DA experiments in Wander et al. *Ecosphere* paper titled "Data assimilation experiments inform monitoring needs for near-term, ecological forecasts in a eutrophic reservoir"

# Instructions to reproduce manuscript + SI figures:

1.  Download or clone github repository to your local computer
2.  Download forecast and score parquet files generated from combined_workflow.R script in the workflows/DA_experiments folder **(LINK TO OTHER ZENODO PUB?**)
3.  Run "BVR_FLARE_ms_figs.R" script in the workflows/DA_experiments folder to reproduce manuscript figures
4.  Run "BVR_FLARE_UC_figs.R" script in the workflows/DA_experiments folder to reproduce Fig. 9 (proportion of IC uncertainty) and SI figures for forecasts run without initial conditions uncertainty
5.  Run "BVR_FLARE_SI_figs.R" script in the workflows/DA_experiments folder to reproduce remaining SI figures
