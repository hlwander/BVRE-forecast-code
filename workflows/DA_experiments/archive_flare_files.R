library(arrow)
library(tidyverse)

site_id_list <- c("bvre")
model_id_list <- c("daily","weekly","fortnightly","monthly")
use_s3 <- FALSE
start_date_list <- seq(as.Date("2020-11-27"), as.Date("2022-02-05"), by="days")

lake_directory <- here::here()

message("Archiving stage 2 NOAA")

s3 <- s3_bucket("drivers/noaa/gefs-v12-reprocess/stage2/parquet/0", endpoint_override = "s3.flare-forecast.org", anonymous = TRUE)

df <- open_dataset(s3, partitioning = c("start_date","site_id")) |> 
  filter(site_id %in% site_id_list, 
         start_date %in% start_date_list)

write_dataset(df, path = file.path(lake_directory, "archive/drivers/drivers/noaa/gefs-v12-reprocess/stage2/parquet/0"), hive_style = FALSE, partitioning = c("start_date","site_id"))

message("Archiving stage 3 NOAA")

s3 <- s3_bucket("drivers/noaa/gefs-v12-reprocess/stage3/parquet", endpoint_override = "s3.flare-forecast.org", anonymous = TRUE)

df <- open_dataset(s3, partitioning = c("site_id")) |> 
  filter(site_id %in% site_id_list)

write_dataset(df, path = file.path(lake_directory, "archive/drivers/drivers/noaa/gefs-v12-reprocess/stage3"), hive_style = FALSE, partitioning = c("site_id"))

setwd(file.path(lake_directory, "archive/drivers"))
files2zip <- fs::dir_ls(recurse = TRUE)
files2zip <- files2zip[stringr::str_detect(files2zip, pattern = "DS_Store", negate = TRUE)][-1]
utils::zip(zipfile = file.path(lake_directory, "archive/drivers"), files = files2zip)

##############

message("Archiving forecast parquets")

if(use_s3){
  s3 <- s3_bucket("forecasts/parquets", endpoint_override = "s3.flare-forecast.org", anonymous = TRUE)
}else{
  s3_all_UC <- file.path(lake_directory, "forecasts/parquets/all_UC")
  s3_IC_off <- file.path(lake_directory, "forecasts/parquets/IC_off")
}

df_all_UC <- open_dataset(s3_all_UC) |> 
  filter(site_id %in% site_id_list,
         model_id %in% model_id_list)

write_dataset(df_all_UC, path = file.path(lake_directory, "archive/forecasts/forecasts/all_UC"), hive_style = TRUE, partitioning = c("site_id","reference_datetime"))

df_IC_off <- open_dataset(s3_IC_off) |> 
  filter(site_id %in% site_id_list,
         model_id %in% model_id_list)

write_dataset(df_IC_off, path = file.path(lake_directory, "archive/forecasts/forecasts/IC_off"), hive_style = TRUE, partitioning = c("site_id","reference_datetime"))

setwd(file.path(lake_directory, "archive/forecasts"))
files2zip <- fs::dir_ls(recurse = TRUE)
files2zip <- files2zip[stringr::str_detect(files2zip, pattern = "DS_Store", negate = TRUE)][-1]
utils::zip(zipfile = file.path(lake_directory, "archive/forecasts"), files = files2zip)

#######


message("Archiving score parquets")

if(use_s3){
  s3 <- s3_bucket("scores/parquets", endpoint_override = "s3.flare-forecast.org", anonymous = TRUE)
}else{
  s3_all_UC <- file.path(lake_directory, "scores/all_UC")
  s3_IC_off <- file.path(lake_directory, "scores/IC_off")
  s3_const_bad <- file.path(lake_directory, "scores/constant_bad_pars")
  s3_tuned_bad <- file.path(lake_directory, "scores/tuned_bad_pars")
}

df_all_UC <- open_dataset(s3_all_UC) |> 
  filter(site_id %in% site_id_list,
         model_id %in% model_id_list)

write_dataset(df_all_UC, path = file.path(lake_directory, "archive/scores/scores/all_UC"), hive_style = TRUE, partitioning = c("site_id","reference_datetime"))

df_IC_off <- open_dataset(s3_IC_off) |> 
  filter(site_id %in% site_id_list,
         model_id %in% model_id_list)

write_dataset(df_IC_off, path = file.path(lake_directory, "archive/scores/scores/IC_off"), hive_style = TRUE, partitioning = c("site_id","reference_datetime"))


df_const_bad <- open_dataset(s3_const_bad) |> 
  filter(site_id %in% site_id_list,
         model_id %in% "daily_no_pars")

write_dataset(df_const_bad, path = file.path(lake_directory, "archive/scores/scores/constant_bad_pars"), hive_style = TRUE, partitioning = c("site_id","reference_datetime"))


df_tuned_bad <- open_dataset(s3_tuned_bad) |> 
  filter(site_id %in% site_id_list,
         model_id %in% "daily_no_pars")

write_dataset(df_tuned_bad, path = file.path(lake_directory, "archive/scores/scores/tuned_bad_pars"), hive_style = TRUE, partitioning = c("site_id","reference_datetime"))


setwd(file.path(lake_directory, "archive/scores"))
files2zip <- fs::dir_ls(recurse = TRUE)
files2zip <- files2zip[stringr::str_detect(files2zip, pattern = "DS_Store", negate = TRUE)][-1]
utils::zip(zipfile = file.path(lake_directory, "archive/scores"), files = files2zip)

#######
message("Archiving targets")

s3 <- file.path(lake_directory, "targets")

fs::dir_create(file.path(lake_directory, "archive/targets"))

file.copy(from = s3, paste0(lake_directory,"/archive/targets"),
          overwrite = TRUE, recursive = TRUE, copy.mode = TRUE)

setwd(file.path(lake_directory, "archive/targets"))
files2zip <- fs::dir_ls(recurse = TRUE)
files2zip <- files2zip[stringr::str_detect(files2zip, pattern = "DS_Store", negate = TRUE)][-1]
utils::zip(zipfile = file.path(lake_directory, "archive/targets"), files = files2zip)
