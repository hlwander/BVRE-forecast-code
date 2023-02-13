pacman::p_load(tidyverse)

noaa_date <- as.character(lubridate::as_date("2021-10-11"))
s3 <- arrow::s3_bucket(bucket = "drivers/noaa/gefs-v12/stage2/parquet/0", endpoint_override = "s3.flare-forecast.org", anonymous = TRUE)
noaa_df <- arrow::open_dataset(s3, partitioning = "reference_date") |> 
  filter(variable == "air_temperature",
         reference_date == noaa_date,
         #horizon <= 16 * 24, 
         site_id %in% c("bvre")) |> 
  select(site_id, reference_datetime, datetime, family, parameter, variable, prediction) |> 
  dplyr::collect() |> 
  mutate(prediction = prediction - 273.15)
noaa_df |> ggplot(aes(x = datetime, y = prediction, group = parameter)) +
  geom_line() +
  labs(y = "NOAA GEFS air temperature") +
  theme_bw()

#calculate vairance for each day

test <- plyr::ddply(noaa_df, c("parameter","datetime"), function(x) {
data.frame(
  mean = mean(x$prediction, na.rm=T)
)
}, .progress = plyr::progress_text(), .parallel = FALSE) 


variance_df <- plyr::ddply(test, c("datetime"), function(x) {
  data.frame(
    sd = sd(x$mean,na.rm=T),
    variance = sd(x$mean, na.rm=T)^2
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 



variance_df |> ggplot(aes(x = datetime, y = variance)) +
  geom_line() +
  labs(y = "NOAA GEFS variance") +
  theme_bw()
