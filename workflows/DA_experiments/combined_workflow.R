library(tidyverse)
library(hydroGOF)
library(lubridate)

set.seed(100)

lake_directory <- here::here()
forecast_site <- "bvre"
configure_run_file <- "configure_run.yml"
config_files <- "configure_flare.yml"
config_set_name <- "DA_experiments"
use_archive <- FALSE

starting_index <- 1

if(use_archive){
  use_s3 <- FALSE
}else{
  Sys.setenv('AWS_DEFAULT_REGION' = 's3',
             'AWS_S3_ENDPOINT' = 'flare-forecast.org',
             'USE_HTTPS' = TRUE)
  use_s3 <- FALSE
}

#DA frequency vectors
daily <- seq.Date(as.Date("2020-11-27"), as.Date("2022-02-01"), by = 1)#changing end to 2022-02-01 because need observations to evaluate forecasts 
date_list <- list(daily = daily,                                      
                  weekly = daily[seq(1, length(daily), 7)],
                  fortnightly = daily[seq(1, length(daily), 14)],
                  monthly = daily[seq(1, length(daily), 30)]) 
models <- names(date_list)

num_forecasts <- 366 #* 3 - 3
#num_forecasts <- 1#19 * 7 + 1
days_between_forecasts <- 1
forecast_horizon <- 35
starting_date <- as_date("2020-11-27") 
second_date <- starting_date + months(1) + days(5) 

start_dates <- as_date(rep(NA, num_forecasts + 1))
end_dates <- as_date(rep(NA, num_forecasts + 1))
start_dates[1] <- starting_date
end_dates[1] <- second_date
for(i in 2:(num_forecasts+1)){
  start_dates[i] <- as_date(end_dates[i-1])
  end_dates[i] <- start_dates[i] + days(days_between_forecasts)
}

j = 1
sites <- "bvre"

sims <- expand.grid(paste0(start_dates,"_",end_dates,"_", forecast_horizon), models)

names(sims) <- c("date","model")

sims$start_dates <- stringr::str_split_fixed(sims$date, "_", 3)[,1]
sims$end_dates <- stringr::str_split_fixed(sims$date, "_", 3)[,2]
sims$horizon <- stringr::str_split_fixed(sims$date, "_", 3)[,3]

sims <- sims |>
  mutate(model = as.character(model)) |>
  select(-date) |>
  distinct_all() |>
  arrange(start_dates)

sims$horizon[1:length(models)] <- 0

message("Generating targets")

#' Source the R files in the repository

source(file.path(lake_directory, "R", "met_qaqc_csv.R"))
source(file.path(lake_directory, "R", "in_situ_qaqc.R"))
source(file.path(lake_directory, "R", "temp_oxy_chla_qaqc.R"))
source(file.path(lake_directory, "R", "extract_CTD.R"))
source(file.path(lake_directory, "R", "extract_secchi.R"))
source(file.path(lake_directory, "R", "extract_nutrients.R"))
source(file.path(lake_directory, "R", "extract_ch4.R"))
source(file.path(lake_directory, "R", "inflow_qaqc.R"))
#source(file.path(lake_directory, "R/old", "TMWB_inflow_model.R"))

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
                                    ctd_fname = file.path(lake_directory, "data_raw", config_obs$ctd_fname),
                                    nutrients_fname =  file.path(lake_directory, "data_raw", config_obs$nutrients_fname),
                                    secchi_fname = file.path(lake_directory, "data_raw", config_obs$secchi_fname),
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



for(i in starting_index:nrow(sims)){
  
  message(paste0("index: ", i))
  message(paste0("     Running model: ", sims$model[i], " "))

  model <- sims$model[i]
  sim_names <- model
  
  config <- FLAREr::set_configuration(configure_run_file,lake_directory, config_set_name = config_set_name, sim_name = sim_names)
  
  cycle <- "00"
  
  if(file.exists(file.path(lake_directory, "restart", sites[j], sim_names, configure_run_file))){
    unlink(file.path(lake_directory, "restart", sites[j], sim_names, configure_run_file))
    if(use_s3){
      FLAREr::delete_restart(site_id = sites[j],
                             sim_name = sim_names,
                             bucket = config$s3$warm_start$bucket,
                             endpoint = config$s3$warm_start$endpoint)
    }
  }
  run_config <- yaml::read_yaml(file.path(lake_directory, "configuration", config_set_name, configure_run_file))
  run_config$configure_flare <- config_files
  run_config$sim_name <- sim_names
  yaml::write_yaml(run_config, file = file.path(lake_directory, "restart", sites[j], sim_names, configure_run_file))
  config <- FLAREr::set_configuration(configure_run_file,lake_directory, config_set_name = config_set_name)
  config$run_config$start_datetime <- as.character(paste0(sims$start_dates[i], " 00:00:00"))
  config$run_config$forecast_start_datetime <- as.character(paste0(sims$end_dates[i], " 00:00:00"))
  config$run_config$forecast_horizon <- sims$horizon[i]
  if(i <= length(models)){
    config$run_config$restart_file <- NA
  }else{
    config$run_config$restart_file <- paste0(config$location$site_id, "-", lubridate::as_date(config$run_config$start_datetime), "-", sim_names, ".nc")
    if(!file.exists(config$run_config$restart_file )){
      warning(paste0("restart file: ", config$run_config$restart_file, " doesn't exist"))
    }
  }
  
  run_config <- config$run_config
  yaml::write_yaml(run_config, file = file.path(lake_directory, "restart", sites[j], sim_names, configure_run_file))
  
  config <- FLAREr::set_configuration(configure_run_file,lake_directory, config_set_name = config_set_name, sim_name = sim_names)
  config$model_settings$model <- model
  config$run_config$sim_name <- sim_names
  config <- FLAREr::get_restart_file(config, lake_directory)
  
  da_freq <- which(names(date_list) == sims$model[i])
  
  met_out <- FLAREr::generate_glm_met_files_arrow(obs_met_file = NULL,
                                                  out_dir = config$file_path$execute_directory,
                                                  forecast_dir = forecast_dir,
                                                  config = config,
                                                  use_s3 = TRUE,
                                                  bucket = config$s3$drivers$bucket,
                                                  endpoint = config$s3$drivers$endpoint)

  met_out$filenames <- met_out$filenames[!stringr::str_detect(met_out$filenames, "31")]
  
  pars_config <- readr::read_csv(file.path(config$file_path$configuration_directory, config$model_settings$par_config_file), col_types = readr::cols())
  obs_config <- readr::read_csv(file.path(config$file_path$configuration_directory, config$model_settings$obs_config_file), col_types = readr::cols())
  states_config <- readr::read_csv(file.path(config$file_path$configuration_directory, config$model_settings$states_config_file), col_types = readr::cols())
  
  #Create observation matrix
  obs <- FLAREr::create_obs_matrix(cleaned_observations_file_long = file.path(config$file_path$qaqc_data_directory, paste0(config$location$site_id, "-targets-insitu.csv")),
                                   obs_config = obs_config,
                                   config)
  
  full_time <- seq(lubridate::as_datetime(config$run_config$start_datetime), lubridate::as_datetime(config$run_config$forecast_start_datetime) + lubridate::days(config$run_config$forecast_horizon), by = "1 day")
  full_time <- as.Date(full_time)
  
  idx <- which(!full_time %in% date_list[[da_freq]])
  obs[1, idx, ] <- NA
  
  #set output directory so each frequency/experiment is saved in a separate folder
  #config$file_path$forecast_output_directory <- file.path(lake_directory,"forecasts","bvre","DA_experiments",names(date_list)[da_freq])
  
  message("Generating forecast")

  #separating this from 2.5 so I can iteratively loop through this scipt for each day I want to generate forecasts
  states_config <- FLAREr::generate_states_to_obs_mapping(states_config, obs_config)
  
  model_sd <- FLAREr::initiate_model_error(config, states_config)
  
  init <- FLAREr::generate_initial_conditions(states_config,
                                              obs_config,
                                              pars_config,
                                              obs,
                                              config,
                                              historical_met_error = met_out$historical_met_error)
  #Run EnKF
  da_forecast_output <- FLAREr::run_da_forecast(states_init = init$states,
                                                pars_init = init$pars,
                                                aux_states_init = init$aux_states_init,
                                                obs = obs,
                                                obs_sd = obs_config$obs_sd,
                                                model_sd = model_sd,
                                                working_directory = config$file_path$execute_directory,
                                                met_file_names = met_out$filenames,
                                                inflow_file_names = NULL,#list.files(inflow_file_dir, pattern='INFLOW-'),
                                                outflow_file_names = NULL,#list.files(inflow_file_dir, pattern='OUTFLOW-'),
                                                config = config,
                                                pars_config = pars_config,
                                                states_config = states_config,
                                                obs_config = obs_config,
                                                management = NULL,
                                                da_method = config$da_setup$da_method,
                                                par_fit_method = config$da_setup$par_fit_method)
  
  saved_file <- FLAREr::write_forecast_netcdf(da_forecast_output = da_forecast_output,
                                              forecast_output_directory = config$file_path$forecast_output_directory,
                                              use_short_filename = TRUE)
  
  forecast_df <- FLAREr::write_forecast_arrow(da_forecast_output = da_forecast_output,
                                              use_s3 = use_s3,
                                              bucket = config$s3$forecasts_parquet$bucket,
                                              endpoint = config$s3$forecasts_parquet$endpoint,
                                              local_directory = file.path(lake_directory, config$s3$forecasts_parquet$bucket))
  
  
  FLAREr::generate_forecast_score_arrow(targets_file = file.path(config$file_path$qaqc_data_directory,paste0(config$location$site_id, "-targets-insitu.csv")),
                                        forecast_df = forecast_df,
                                        use_s3 = use_s3,
                                        bucket = config$s3$scores$bucket,
                                        endpoint = config$s3$scores$endpoint,
                                        local_directory = file.path(lake_directory, config$s3$scores$bucket),
                                        variable_types = c("state","parameter"))
  
  #rm(da_forecast_output)
  #gc()
  message("Generating plot")
  FLAREr::plotting_general_2(file_name = saved_file,
                             target_file = file.path(config$file_path$qaqc_data_directory, paste0(config$location$site_id, "-targets-insitu.csv")),
                             ncore = 2,
                             obs_csv = FALSE)
  
  FLAREr::put_forecast(saved_file, eml_file_name, config)
  
  #unlink(saved_file)
  
  #unlink(config$run_config$restart_file)
  
  rm(da_forecast_output)
  gc()
  
  sink('last_completed_index.txt')
  print(i)
  sink()
  
}
