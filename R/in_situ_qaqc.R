in_situ_qaqc <- function(insitu_obs_fname,
                         data_location,
                         maintenance_file,
                         ctd_fname,
                         nutrients_fname,
                         secchi_fname,
                         cleaned_insitu_file,
                         site_id,
                         config){

  print("QAQC BVR sensors")

  d <- temp_oxy_chla_qaqc(realtime_file = insitu_obs_fname[1],
                          qaqc_file = insitu_obs_fname[2],
                          offset_file = insitu_obs_fname[3],
                          maintenance_file = maintenance_file,
                          config = config)

  if(exists("ctd_fname")){
    if(!is.na(ctd_fname)){
      print("QAQC CTD")
      d_ctd <- extract_CTD(fname = file.path(config$file_path$data_directory,config$ctd_fname),
                           input_file_tz = "EST",
                           focal_depths = config$focal_depths,
                           config = config)
      d <- rbind(d,d_ctd)
    }
  }


  if(exists("nutrients_fname")){
    if(!is.na(nutrients_fname)){
      print("QAQC Nutrients")
      d_nutrients <- extract_nutrients(fname = file.path(config$file_path$data_directory,config$nutrients_fname),
                                       input_file_tz = "EST",
                                       focal_depths = config$focal_depths)
      d <- rbind(d,d_nutrients)
    }
  }


  if(exists("ch4_fname")){
    if(!is.na(ch4_fname)){
      print("QAQC CH4")
      d_ch4 <- extract_ch4(fname = file.path(config$file_path$data_directory,config$ch4_fname),
                           input_file_tz = "EST",
                           focal_depths = config$focal_depths)
      d <- rbind(d,d_ch4)
    }
  }
  
  #drop NA rows
  d <- d[!is.na(d$observed),]
  
  first_day  <- lubridate::as_datetime(paste0(lubridate::as_date(min(d$time)), " ", config$averaging_period_starting_hour))
  first_day <- lubridate::force_tz(first_day, tzone = "UTC")

  last_day <- lubridate::as_datetime(paste0(lubridate::as_date(max(d$time)), " ", config$averaging_period_starting_hour))
  last_day <- lubridate::force_tz(last_day, tzone = "UTC")

  full_time_local <- seq(first_day, last_day, by = "1 day")

  d_clean <- NULL


  for(i in 1:length(config$target_variable)){
    print(paste0("Extracting ",config$target_variable[i]))
    #depth_breaks <- sort(c(bins1, bins2))
    time_breaks <- seq(first_day, last_day, by = config$averaging_period[i])

    d_curr <- d %>%
      dplyr::filter(variable == config$target_variable[i],
                    method %in% config$measurement_methods[[i]]) %>% #NOTE, ADDED CTD AS METHOD SO WE HAVE TEMP AND DO DATA BEFORE 2020
      dplyr::mutate(time_class = cut(time, breaks = time_breaks, labels = FALSE)) %>%
      dplyr::group_by(time_class, depth) %>%
      dplyr::summarize(observed = mean(observed, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(time = time_breaks[time_class]) %>%
      dplyr::mutate(variable = config$target_variable[i]) %>%
      dplyr::select(time, depth, observed, variable) %>%
      dplyr::mutate(date = lubridate::as_date(time))

    if(config$averaging_period[i] == "1 hour"){
      d_curr <- d_curr %>%
      dplyr::mutate(hour = lubridate::hour(time)) %>%
      dplyr::filter(hour == lubridate::hour(first_day))  #getting rid of this filter because I need more observations!
    }else{
      d_curr <- d_curr %>%
        dplyr::mutate(hour = NA) %>%
        dplyr::mutate(hour = as.numeric(hour))
    }

    d_curr <- d_curr %>% dplyr::select(-date)

    d_clean <- rbind(d_clean,  d_curr) 
  }

  d_clean <- d_clean %>% tidyr::drop_na(observed)

  if(!is.na(secchi_fname)){

    d_secchi <- extract_secchi(fname = secchi_fname,
                               input_file_tz = "EST",
                               focal_depths = config$focal_depths)

    d_secchi <- d_secchi %>%
      dplyr::mutate(time = lubridate::as_date(time)) %>%
      dplyr::mutate(hour = 0) 

    d_clean <- rbind(d_clean,d_secchi) 
  }

  d_clean$site_id <- "bvre"
  
  d_clean <- d_clean %>% 
    dplyr::select(time, site_id, depth, observed, variable) |> 
    dplyr::rename(datetime = time,
                  observation = observed)
  
  d_clean$observation <- round(d_clean$observation, digits = 4)
  
  if(!dir.exists(dirname(cleaned_insitu_file))){
    dir.create(dirname(cleaned_insitu_file), recursive = TRUE)
  }
  
  readr::write_csv(d_clean, cleaned_insitu_file)
  
  return(cleaned_insitu_file)

}
       
