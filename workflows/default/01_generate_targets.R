library(tidyverse)
library(lubridate)

message("Beginning generate targets")

#' Set the lake directory to the repository directory

lake_directory <- here::here()
config_set_name <- "default"

Sys.setenv("AWS_DEFAULT_REGION" = "s3",
           "AWS_S3_ENDPOINT" = "flare-forecast.org")

#' Source the R files in the repository
files.sources <- list.files(file.path(lake_directory, "R"), full.names = TRUE)
sapply(files.sources[grepl(".R$", files.sources)], source)

#' Generate the `config_obs` object and create directories if necessary

config_obs <- FLAREr::initialize_obs_processing(lake_directory, observation_yml = "observation_processing.yml", config_set_name = config_set_name)
dir.create(file.path(lake_directory, "targets", config_obs$site_id), showWarnings = FALSE)
#' Clone or pull from data repositories

FLAREr::get_git_repo(lake_directory,
                     directory = config_obs$realtime_insitu_location,
                     git_repo = "https://github.com/FLARE-forecast/BVRE-data.git")

#' Download files from EDI

FLAREr::get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/725/2/026a6e2cca8bdf18720d6a10d8860e3d",
                     file = config_obs$insitu_obs_fname[2],
                     lake_directory)

FLAREr::get_edi_file(edi_https = "https://pasta.lternet.edu/package/data/eml/edi/725/2/8c0d1d8ea078d274c252cd362a500d26",
                     file = config_obs$insitu_obs_fname[3],
                     lake_directory)


#' Clean up observed insitu measurements
cleaned_insitu_file <- in_situ_qaqc(insitu_obs_fname = file.path(lake_directory,"data_raw", config_obs$insitu_obs_fname),
                                    data_location = file.path(lake_directory,"data_raw"),
                                    maintenance_file = file.path(lake_directory, "data_raw", config_obs$maintenance_file),
                                    ctd_fname = NA,
                                    nutrients_fname =  NA,
                                    secchi_fname = NA,
                                    cleaned_insitu_file = file.path(lake_directory,"targets", config_obs$site_id, paste0(config_obs$site_id,"-targets-insitu.csv")),
                                    site_id = config_obs$site_id,
                                    config = config_obs)

#' Move targets to s3 bucket

message("Successfully generated targets")

FLAREr::put_targets(site_id = config_obs$site_id,
                    cleaned_insitu_file,
                    cleaned_met_file,
                    use_s3 = FALSE)

message("Successfully moved targets to s3 bucket")
