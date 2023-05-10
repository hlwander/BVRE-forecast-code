library(arrow)
library(tidyverse)

site_id_list <- c("bvre")
model_id_list <- c("Daily","Weekly","Fortnightly","Monthly")
use_s3 <- FALSE
start_date_list <- seq(as.Date("2020-11-27"), as.Date("2022-02-05"), by="days")

lake_directory <- here::here()

fs::dir_create(file.path(lake_directory, "archive/forecasts/forecasts"))
fs::dir_create(file.path(lake_directory, "archive/drivers/drivers"))
fs::dir_create(file.path(lake_directory, "archive/scores/scores"))
fs::dir_create(file.path(lake_directory, "archive/targets/targets"))

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
utils::zip(zipfile = file.path(lake_directory, "archive/drivers"), files = files2zip)

##############

message("Archiving forecast parquets")

if(use_s3){
  s3 <- s3_bucket("forecasts/parquets", endpoint_override = "s3.flare-forecast.org", anonymous = TRUE)
}else{
  s3 <- file.path(lake_directory, "forecasts/all_UC")
}

df <- open_dataset(s3) |> 
  filter(site_id %in% site_id_list,
         model_id %in% model_id_list)

write_dataset(df, path = file.path(lake_directory, "archive/forecasts"), hive_style = TRUE, partitioning = c("site_id","reference_datetime"))

setwd(file.path(lake_directory, "archive/forecasts"))
files2zip <- fs::dir_ls(recurse = TRUE)
utils::zip(zipfile = file.path(lake_directory, "archive/forecasts"), files = files2zip)

#######


message("Archiving score parquets")

if(use_s3){
  s3 <- s3_bucket("scores/parquets", endpoint_override = "s3.flare-forecast.org", anonymous = TRUE)
}else{
  s3 <- file.path(lake_directory, "scores")
}

df <- open_dataset(s3) |> 
  filter(site_id %in% site_id_list,
         model_id %in% model_id_list)

write_dataset(df, path = file.path(lake_directory, "archive/scores/scores"), hive_style = TRUE, partitioning = c("site_id","reference_datetime"))

setwd(file.path(lake_directory, "archive/scores"))
files2zip <- fs::dir_ls(recurse = TRUE)
utils::zip(zipfile = file.path(lake_directory, "archive/scores"), files = files2zip)

#######
message("Archiving targets")

s3 <- file.path(lake_directory, "targets")

file.copy(from = s3, paste0(lake_directory,"/archive/targets/targets/"),
          overwrite = TRUE, recursive = TRUE, copy.mode = TRUE)

setwd(file.path(lake_directory, "archive/targets"))
files2zip <- fs::dir_ls(recurse = TRUE)
utils::zip(zipfile = file.path(lake_directory, "archive/targets"), files = files2zip)
