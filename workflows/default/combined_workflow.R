library(tidyverse)
library(lubridate)
lake_directory <- here::here()
setwd(lake_directory)
forecast_site <<- "bvre"
configure_run_file <<- "configure_run.yml"
update_run_config <<- TRUE
config_set_name <<- "default"

FLAREr::ignore_sigpipe()

noaa_ready <- TRUE

while(noaa_ready){
    
  message("Generating targets")
  source(file.path("workflows", config_set_name, "01_generate_targets.R"))
  
  setwd(lake_directory)
  
  message("Generating inflow forecast")
  source(file.path("workflows", config_set_name, "02_run_inflow_forecast.R"))
  
  setwd(lake_directory)
  
  message("Generating forecast")
  source(file.path("workflows", config_set_name, "03_run_flarer_forecast.R"))
  
  setwd(lake_directory)
    
  RCurl::url.exists("https://hc-ping.com/8b5c849d-a5a6-4d44-980c-676472cf3c70", timeout = 5)
  
  noaa_ready <- FLAREr::check_noaa_present_arrow(lake_directory,
                                               configure_run_file,
                                               config_set_name = config_set_name)
}
