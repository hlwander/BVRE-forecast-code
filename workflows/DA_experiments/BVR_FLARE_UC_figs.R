#UC analysis figs
#29 Dec 2022

#load libraries
pacman::p_load(ggplot2,tidyverse,FSA,rcompanion)

# filtering function - turns outliers into NAs to be removed
filter_lims <- function(x){
  l <- boxplot.stats(x)$stats[1]
  u <- boxplot.stats(x)$stats[5]
  
  for (i in 1:length(x)){
    x[i] <- ifelse(x[i]>l & x[i]<u, x[i], NA)
  }
  return(x)
}


#set wd
lake_directory <- here::here()
setwd(lake_directory)

#read in all forecasts 
score_dir <- arrow::SubTreeFileSystem$create(file.path(lake_directory,"scores/UC"))
all_DA_forecasts <- arrow::open_dataset(score_dir) |> collect() |>   
  filter(!is.na(observation), variable == "temperature",horizon > 0)

#need to round horizon becuase they are in decimal form for some reason...
all_DA_forecasts$horizon <- ceiling(all_DA_forecasts$horizon)

#round depths up to nearest m 
all_DA_forecasts$depth <- ceiling(all_DA_forecasts$depth)

#add a group number so that I can average horizons later on
all_DA_forecasts <- all_DA_forecasts %>% 
  mutate(group = case_when(all_DA_forecasts$horizon <= 5 ~ "1-5",
                           all_DA_forecasts$horizon <=10 & all_DA_forecasts$horizon > 5 ~ "6-10",
                           all_DA_forecasts$horizon <=15 & all_DA_forecasts$horizon > 10 ~ "11-15",
                           all_DA_forecasts$horizon <=20 & all_DA_forecasts$horizon > 15 ~ "16-20",
                           all_DA_forecasts$horizon <=25 & all_DA_forecasts$horizon > 20 ~ "21-25",
                           all_DA_forecasts$horizon <=30 & all_DA_forecasts$horizon > 25 ~ "26-30",
                           all_DA_forecasts$horizon <=36 & all_DA_forecasts$horizon > 30 ~ "31-35"))

strat_date<- "2021-11-07"

#add stratified vs mixed col
all_DA_forecasts$phen <- ifelse(all_DA_forecasts$datetime <= as.POSIXct(strat_date) & 
                                  all_DA_forecasts$datetime >="2021-03-13","Stratified", "Mixed")

#remove n=6 days with ice-cover 
all_DA_forecasts <- all_DA_forecasts[!(as.Date(all_DA_forecasts$datetime) %in% c(as.Date("2021-01-10"), as.Date("2021-01-11"),as.Date("2021-01-30"),
                                                                                 as.Date("2021-02-13"),as.Date("2021-02-14"),as.Date("2021-02-15"))),]

#change model_id to be all uppercase
all_DA_forecasts$model_id <- str_to_title(all_DA_forecasts$model_id)

#only keep 2021 data
all_DA_forecasts <- all_DA_forecasts[all_DA_forecasts$datetime<="2021-12-31",]

#------------------------------------------------------------------------------#
#calculate forecast skill metrics

#forecast skill for each depth and horizon
forecast_skill_depth_horizon <-  plyr::ddply(all_DA_forecasts, c("depth","horizon","phen", "model_id"), function(x) { #phen instead of datetime?
  data.frame(
    RMSE = sqrt(mean((x$mean - x$observation)^2, na.rm = TRUE)),
    MAE = mean(abs(x$mean - x$observation), na.rm = TRUE),
    pbias = 100 * (sum(x$mean - x$observation, na.rm = TRUE) / sum(x$observation, na.rm = TRUE)),
    CRPS = verification::crps(x$observation, as.matrix(x[, c(7,9)]))$CRPS,
    variance = (mean(x$sd))^2
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 


#order DA frequencies
forecast_skill_depth_horizon$model_id <- factor(forecast_skill_depth_horizon$model_id, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))


#df with averaged forecast skill for all days (group by horizon, DA, and phen)
forecast_horizon_avg <- plyr::ddply(all_DA_forecasts, c("horizon", "model_id", "phen"), function(x) {
  data.frame(
    RMSE = sqrt(mean((x$mean - x$observation)^2, na.rm = TRUE)),
    MAE = mean(abs(x$mean - x$observation), na.rm = TRUE),
    pbias = 100 * (sum(x$mean - x$observation, na.rm = TRUE) / sum(x$observation, na.rm = TRUE)),
    CRPS = verification::crps(x$observation, as.matrix(x[, c(7,9)]))$CRPS
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

#order DA frequencies
forecast_horizon_avg$model_id <- factor(forecast_horizon_avg$model_id, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))

#df with averaged forecast skill for all days (group by depth, DA, and phen)
forecast_depth_avg <- plyr::ddply(all_DA_forecasts, c("depth", "model_id", "phen"), function(x) {
  data.frame(
    RMSE = sqrt(mean((x$mean - x$observation)^2, na.rm = TRUE)),
    MAE = mean(abs(x$mean - x$observation), na.rm = TRUE),
    pbias = 100 * (sum(x$mean - x$observation, na.rm = TRUE) / sum(x$observation, na.rm = TRUE)),
    CRPS = verification::crps(x$observation, as.matrix(x[, c(7,9)]))$CRPS,
    variance = (mean(x$sd))^2
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

#order DA frequencies
forecast_depth_avg$model_id <- factor(forecast_depth_avg$model_id, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))


#------------------------------------------------------------------------------#
#FIGURES
cb_friendly_2 <- c("#8C510A", "#BF812D", "#C7EAE5", "#35978F") #"#DFC27D", "#DEDEDE",

#FIGURE 5: DA frequency boxplots for 1, 5, and 9m and 1, 7, and 35-day horizons

#creating smaller dataset for kw test w/ 1,5,9m and 1,7,35 days
kw_horizons <- forecast_skill_depth_horizon[forecast_skill_depth_horizon$depth %in% c(1,5,9) & forecast_skill_depth_horizon$horizon %in% c(1,7,35),]

#kruskal wallis and dunn tests for 1m stratified rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==1])
dunn_strat_1m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1] ~ 
                            kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==1])
rslt_strat_1m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_1m$res, threshold = 0.05)$Letter[c(1,4,2,3)])
median_strat_1m_daily <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Daily"]),
                           median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Daily"]),
                           median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Daily"]))
median_strat_1m_weekly <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Weekly"]),
                            median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Weekly"]),
                            median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Weekly"]))
median_strat_1m_fortnightly <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Fortnightly"]),
                                 median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Fortnightly"]),
                                 median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Fortnightly"]))
median_strat_1m_monthly <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Monthly"]),
                             median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Monthly"]),
                             median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Monthly"]))

#kruskal wallis and dunn tests for 5m stratified rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==5]))
dunn_strat_5m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5] ~ 
                            as.factor(kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==5]))
rslt_strat_5m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_5m$res, threshold = 0.05)$Letter[c(1,4,2,3)])
median_strat_5m_daily <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Daily"]),
                           median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Daily"]),
                           median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Daily"]))
median_strat_5m_weekly <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Weekly"]),
                            median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Weekly"]),
                            median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Weekly"]))
median_strat_5m_fortnightly <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Fortnightly"]),
                                 median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Fortnightly"]),
                                 median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Fortnightly"]))
median_strat_5m_monthly <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Monthly"]),
                             median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Monthly"]),
                             median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Monthly"]))

#kruskal wallis and dunn tests for 9m stratified rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==9]))
dunn_strat_9m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9] ~ 
                            as.factor(kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==9]))
rslt_strat_9m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_9m$res, threshold = 0.05)$Letter[c(1,4,2,3)])
median_strat_9m_daily <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Daily"]),
                           median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Daily"]),
                           median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Daily"]))
median_strat_9m_weekly <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Weekly"]),
                            median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Weekly"]),
                            median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Weekly"]))
median_strat_9m_fortnightly <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Fortnightly"]),
                                 median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Fortnightly"]),
                                 median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Fortnightly"]))
median_strat_9m_monthly <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Monthly"]),
                             median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Monthly"]),
                             median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Monthly"]))

#kruskal wallis and dunn tests for 1m mixed rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==1]))
dunn_mix_1m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1] ~ 
                          as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==1]))
rslt_mix_1m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_1m$res, threshold = 0.05)$Letter[c(1,4,2,3)])
median_mix_1m_daily <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Daily"]),
                         median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Daily"]),
                         median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Daily"]))
median_mix_1m_weekly <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Weekly"]),
                          median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Weekly"]),
                          median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Weekly"]))
median_mix_1m_fortnightly <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Fortnightly"]),
                               median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Fortnightly"]),
                               median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Fortnightly"]))
median_mix_1m_monthly <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Monthly"]),
                           median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Monthly"]),
                           median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Monthly"]))

#kruskal wallis and dunn tests for 5m mixed rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==5]))
dunn_mix_5m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5] ~ 
                          as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==5]))
rslt_mix_5m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_5m$res, threshold = 0.05)$Letter[c(1,4,2,3)])
median_mix_5m_daily <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Daily"]),
                         median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Daily"]),
                         median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Daily"]))
median_mix_5m_weekly <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Weekly"]),
                          median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Weekly"]),
                          median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Weekly"]))
median_mix_5m_fortnightly <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Fortnightly"]),
                               median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Fortnightly"]),
                               median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Fortnightly"]))
median_mix_5m_monthly <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Monthly"]),
                           median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Monthly"]),
                           median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Monthly"]))

#kruskal wallis and dunn tests for 9m mixed rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==9]))
dunn_mix_9m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9] ~ 
                          as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==9]))
rslt_mix_9m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_9m$res, threshold = 0.05)$Letter[c(1,4,2,3)])
median_mix_9m_daily <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Daily"]),
                         median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Daily"]),
                         median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Daily"]))
median_mix_9m_weekly <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Weekly"]),
                          median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Weekly"]),
                          median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Weekly"]))
median_mix_9m_fortnightly <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Fortnightly"]),
                               median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Fortnightly"]),
                               median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Fortnightly"]))
median_mix_9m_monthly <- c(median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Monthly"]),
                           median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Monthly"]),
                           median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Monthly"]))

#create table to export mixed and stratified p-vals across depths (aggregated over horizon)
pvals_horizon_aggregated <- data.frame("Depth_m"= c(rep(1,6),rep(5,6),rep(9,6)), 
                                       "Comparison" = c(dunn_mix_1m$res$Comparison,dunn_mix_5m$res$Comparison,dunn_mix_9m$res$Comparison),
                                       #"Z" = c(dunn_mix_1m$res$Z,dunn_mix_5m$res$Z,dunn_mix_9m$res$Z),
                                       "Mixed_pvalue" = c(dunn_mix_1m$res$P.adj,dunn_mix_5m$res$P.adj,dunn_mix_9m$res$P.adj),
                                       "Stratified_pvalue" = c(dunn_strat_1m$res$P.adj,dunn_strat_5m$res$P.adj,dunn_strat_9m$res$P.adj))

#write.csv(pvals_horizon_aggregated,file.path(lake_directory,"analysis/data/UC_pvals_horizon_aggregatred.csv"),row.names = FALSE)

#median RMSE table
median_RMSE_horizon <- data.frame("Depth_m" = c(rep(1,24),rep(5,24),rep(9,24)),
                                  "model_id" = rep(c(rep("Daily",3), rep("Weekly",3),rep("Fortnightly",3),rep("Monthly",3)),6),
                                  "Horizon_days" = rep(c(1,7,35),24),
                                  "TempDynamics" = rep(c(rep("Mixed",12),rep("Stratified",12)),3),
                                  "RMSE_C" = c(median_mix_1m_daily,median_mix_1m_weekly,median_mix_1m_fortnightly,median_mix_1m_monthly,
                                               median_strat_1m_daily,median_strat_1m_weekly,median_strat_1m_fortnightly,median_strat_1m_monthly,
                                               median_mix_5m_daily,median_mix_5m_weekly,median_mix_5m_fortnightly,median_mix_5m_monthly,
                                               median_strat_5m_daily,median_strat_5m_weekly,median_strat_5m_fortnightly,median_strat_5m_monthly,
                                               median_mix_9m_daily,median_mix_9m_weekly,median_mix_9m_fortnightly,median_mix_9m_monthly,
                                               median_strat_9m_daily,median_strat_9m_weekly,median_strat_9m_fortnightly,median_strat_9m_monthly))
#write.csv(median_RMSE_horizon,file.path(lake_directory,"analysis/data/UC_median_RMSE_depth_horizon_DA.csv"),row.names = FALSE)

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$TempDynamics=="Mixed"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$TempDynamics=="Stratified"])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==1])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==7])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==35])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==1])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==5])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==9])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==1 & median_RMSE_horizon$TempDynamics=="Mixed"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==5 & median_RMSE_horizon$TempDynamics=="Mixed"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==9 & median_RMSE_horizon$TempDynamics=="Mixed"])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==1 & median_RMSE_horizon$TempDynamics=="Stratified"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==5 & median_RMSE_horizon$TempDynamics=="Stratified"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==9 & median_RMSE_horizon$TempDynamics=="Stratified"])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==1 & median_RMSE_horizon$TempDynamics=="Mixed"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==7 & median_RMSE_horizon$TempDynamics=="Mixed"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==35 & median_RMSE_horizon$TempDynamics=="Mixed"])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==1 & median_RMSE_horizon$Depth_m==1 & median_RMSE_horizon$TempDynamics=="Mixed"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==1 & median_RMSE_horizon$Depth_m==5 & median_RMSE_horizon$TempDynamics=="Mixed"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==1 & median_RMSE_horizon$Depth_m==9 & median_RMSE_horizon$TempDynamics=="Mixed"])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==35 & median_RMSE_horizon$Depth_m==1 & median_RMSE_horizon$TempDynamics=="Mixed"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==35 & median_RMSE_horizon$Depth_m==5 & median_RMSE_horizon$TempDynamics=="Mixed"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==35 & median_RMSE_horizon$Depth_m==9 & median_RMSE_horizon$TempDynamics=="Mixed"])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==1 & median_RMSE_horizon$Depth_m==1 & median_RMSE_horizon$TempDynamics=="Stratified"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==1 & median_RMSE_horizon$Depth_m==5 & median_RMSE_horizon$TempDynamics=="Stratified"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==1 & median_RMSE_horizon$Depth_m==9 & median_RMSE_horizon$TempDynamics=="Stratified"])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==7 & median_RMSE_horizon$Depth_m==1 & median_RMSE_horizon$TempDynamics=="Stratified"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==7 & median_RMSE_horizon$Depth_m==5 & median_RMSE_horizon$TempDynamics=="Stratified"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==7 & median_RMSE_horizon$Depth_m==9 & median_RMSE_horizon$TempDynamics=="Stratified"])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==35 & median_RMSE_horizon$Depth_m==1 & median_RMSE_horizon$TempDynamics=="Stratified"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==35 & median_RMSE_horizon$Depth_m==5 & median_RMSE_horizon$TempDynamics=="Stratified"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==35 & median_RMSE_horizon$Depth_m==9 & median_RMSE_horizon$TempDynamics=="Stratified"])


mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==1 & median_RMSE_horizon$TempDynamics=="Stratified"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==7 & median_RMSE_horizon$TempDynamics=="Stratified"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==35 & median_RMSE_horizon$TempDynamics=="Stratified"])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==1 & median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Daily"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==1 & median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Weekly"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==1 & median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Fortnightly"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==1 & median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Monthly"])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==5 & median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Daily"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==5 & median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Weekly"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==5 & median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Fortnightly"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==5 & median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Monthly"])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==9 & median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Daily"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==9 & median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Weekly"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==9 & median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Fortnightly"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==9 & median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Monthly"])


#calcualte percent difference of 35 and 1 day RMSE for mixed vs. stratified
(mean(median_RMSE_horizon$RMSE[median_RMSE_horizon$Horizon_days==35 & median_RMSE_horizon$TempDynamics=="Mixed"]) -
    mean(median_RMSE_horizon$RMSE[median_RMSE_horizon$Horizon_days==1 & median_RMSE_horizon$TempDynamics=="Mixed"])) /
  mean(median_RMSE_horizon$RMSE[median_RMSE_horizon$TempDynamics=="Mixed"]) *100

(mean(median_RMSE_horizon$RMSE[median_RMSE_horizon$Horizon_days==35 & median_RMSE_horizon$TempDynamics=="Stratified"]) -
    mean(median_RMSE_horizon$RMSE[median_RMSE_horizon$Horizon_days==1 & median_RMSE_horizon$TempDynamics=="Stratified"])) /
  mean(median_RMSE_horizon$RMSE[median_RMSE_horizon$TempDynamics=="Stratified"]) *100

#mixed vs stratified rmse across depths
(mean(median_RMSE_horizon$RMSE[median_RMSE_horizon$Depth_m==1 & median_RMSE_horizon$TempDynamics=="Stratified"]) - 
    mean(median_RMSE_horizon$RMSE[median_RMSE_horizon$Depth_m==1 & median_RMSE_horizon$TempDynamics=="Mixed"])) /
  mean(median_RMSE_horizon$RMSE[median_RMSE_horizon$Depth_m==1]) *100

(mean(median_RMSE_horizon$RMSE[median_RMSE_horizon$Depth_m==5 & median_RMSE_horizon$TempDynamics=="Stratified"]) - 
    mean(median_RMSE_horizon$RMSE[median_RMSE_horizon$Depth_m==5 & median_RMSE_horizon$TempDynamics=="Mixed"])) /
  mean(median_RMSE_horizon$RMSE[median_RMSE_horizon$Depth_m==5]) *100

(mean(median_RMSE_horizon$RMSE[median_RMSE_horizon$Depth_m==9 & median_RMSE_horizon$TempDynamics=="Stratified"]) -
    mean(median_RMSE_horizon$RMSE[median_RMSE_horizon$Depth_m==9 & median_RMSE_horizon$TempDynamics=="Mixed"])) /
  mean(median_RMSE_horizon$RMSE[median_RMSE_horizon$Depth_m==9]) *100


#kruskal wallis and dunn tests for 1m stratified rmse across horizon and DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Daily"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Daily"]))
dunn_strat_daily_1m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Daily"] ~ 
                                  as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Daily"]))
rslt_strat_daily_1m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_daily_1m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_strat_daily_1m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==1], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==7], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==35], probs=.75))


kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly"]))
dunn_strat_weekly_1m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly"] ~ 
                                   as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly"]))
rslt_strat_weekly_1m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_weekly_1m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_strat_weekly_1m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==1], probs=.75),
                              quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==7], probs=.75),
                              quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly"]))
dunn_strat_fortnightly_1m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly"] ~ 
                                        as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly"]))
rslt_strat_fortnightly_1m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_fortnightly_1m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_strat_fortnightly_1m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==1], probs=.75),
                                   quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==7], probs=.75),
                                   quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly"]))
dunn_strat_monthly_1m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly"] ~ 
                                    as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly"]))
rslt_strat_monthly_1m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_monthly_1m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_strat_monthly_1m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly" & kw_horizons$horizon ==1], probs=.75),
                               quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly" & kw_horizons$horizon ==7], probs=.75),
                               quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly" & kw_horizons$horizon ==35], probs=.75))


#kruskal wallis and dunn tests for 1-day ahead stratified rmse across depth and DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Daily"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Daily"]))
dunn_strat_daily_1day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Daily"] ~ 
                                    as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Daily"]))
rslt_strat_daily_1day=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_daily_1day$res, threshold = 0.05)$Letter)
rslt_strat_daily_1day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Daily" & kw_horizons$depth ==1], probs=.75),
                               quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Daily" & kw_horizons$depth ==5], probs=.75),
                               quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Daily" & kw_horizons$depth ==9], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Weekly"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Weekly"]))
dunn_strat_weekly_1day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Weekly"] ~ 
                                     as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Weekly"]))
rslt_strat_weekly_1day=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_weekly_1day$res, threshold = 0.05)$Letter)
rslt_strat_weekly_1day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Weekly" & kw_horizons$depth ==1], probs=.75),
                                quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Weekly" & kw_horizons$depth ==5], probs=.75),
                                quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Weekly" & kw_horizons$depth ==9], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Fortnightly"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Fortnightly"]))
dunn_strat_fortnightly_1day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Fortnightly"] ~ 
                                          as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Fortnightly"]))
rslt_strat_fortnightly_1day=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_fortnightly_1day$res, threshold = 0.05)$Letter)
rslt_strat_fortnightly_1day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Fortnightly" & kw_horizons$depth ==1], probs=.75),
                                     quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Fortnightly" & kw_horizons$depth ==5], probs=.75),
                                     quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Fortnightly" & kw_horizons$depth ==9], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Monthly"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Monthly"]))
dunn_strat_monthly_1day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Monthly"] ~ 
                                      as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Monthly"]))
rslt_strat_monthly_1day=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_monthly_1day$res, threshold = 0.05)$Letter)
rslt_strat_monthly_1day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Monthly" & kw_horizons$depth ==1], probs=.75),
                                 quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Monthly" & kw_horizons$depth ==5], probs=.75),
                                 quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Monthly" & kw_horizons$depth ==9], probs=.75))


#kruskal wallis and dunn tests for 5m stratified rmse across horizon and DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Daily"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Daily"]))
dunn_strat_daily_5m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Daily"] ~ 
                                  as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Daily"]))
rslt_strat_daily_5m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_daily_5m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_strat_daily_5m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==1], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==7], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly"]))
dunn_strat_weekly_5m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly"] ~ 
                                   as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly"]))
rslt_strat_weekly_5m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_weekly_5m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_strat_weekly_5m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==1], probs=.75),
                              quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==7], probs=.75),
                              quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly"]))
dunn_strat_fortnightly_5m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly"] ~ 
                                        as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly"]))
rslt_strat_fortnightly_5m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_fortnightly_5m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_strat_fortnightly_5m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==1], probs=.75),
                                   quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==7], probs=.75),
                                   quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Monthly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Monthly"]))
dunn_strat_monthly_5m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Monthly"] ~ 
                                    as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Monthly"]))
rslt_strat_monthly_5m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_monthly_5m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_strat_monthly_5m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Monthly" & kw_horizons$horizon ==1], probs=.75),
                               quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Monthly" & kw_horizons$horizon ==7], probs=.75),
                               quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Monthly" & kw_horizons$horizon ==35], probs=.75))

#kruskal wallis and dunn tests for 7day ahead stratified rmse across depths and DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Daily"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Daily"]))
dunn_strat_daily_7day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Daily"] ~ 
                                    as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Daily"]))
rslt_strat_daily_7day=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_daily_7day$res, threshold = 0.05)$Letter)
rslt_strat_daily_7day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Daily" & kw_horizons$depth ==1], probs=.75),
                               quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Daily" & kw_horizons$depth ==5], probs=.75),
                               quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Daily" & kw_horizons$depth ==9], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Weekly"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Weekly"]))
dunn_strat_weekly_7day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Weekly"] ~ 
                                     as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Weekly"]))
rslt_strat_weekly_7day=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_weekly_7day$res, threshold = 0.05)$Letter)
rslt_strat_weekly_7day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Weekly" & kw_horizons$depth ==1], probs=.75),
                                quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Weekly" & kw_horizons$depth ==5], probs=.75),
                                quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Weekly" & kw_horizons$depth ==9], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Fortnightly"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Fortnightly"]))
dunn_strat_fortnightly_7day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Fortnightly"] ~ 
                                          as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Fortnightly"]))
rslt_strat_fortnightly_7day=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_fortnightly_7day$res, threshold = 0.05)$Letter)
rslt_strat_fortnightly_7day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Fortnightly" & kw_horizons$depth ==1], probs=.75),
                                     quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Fortnightly" & kw_horizons$depth ==5], probs=.75),
                                     quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Fortnightly" & kw_horizons$depth ==9], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Monthly"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Monthly"]))
dunn_strat_monthly_7day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Monthly"] ~ 
                                      as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Monthly"]))
rslt_strat_monthly_7day=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_monthly_7day$res, threshold = 0.05)$Letter)
rslt_strat_monthly_7day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Monthly" & kw_horizons$depth ==1], probs=.75),
                                 quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Monthly" & kw_horizons$depth ==5], probs=.75),
                                 quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7 & kw_horizons$model_id == "Monthly" & kw_horizons$depth ==9], probs=.75))


#kruskal wallis and dunn tests for 9m stratified rmse across horizon and DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Daily"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Daily"]))
dunn_strat_daily_9m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Daily"] ~ 
                                  as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Daily"]))
rslt_strat_daily_9m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_daily_9m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_strat_daily_9m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==1], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==7], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly"]))
dunn_strat_weekly_9m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly"] ~ 
                                   as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly"]))
rslt_strat_weekly_9m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_weekly_9m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_strat_weekly_9m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==1], probs=.75),
                              quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==7], probs=.75),
                              quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly"]))
dunn_strat_fortnightly_9m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly"] ~ 
                                        as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly"]))
rslt_strat_fortnightly_9m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_fortnightly_9m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_strat_fortnightly_9m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==1], probs=.75),
                                   quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==7], probs=.75),
                                   quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Monthly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Monthly"]))
dunn_strat_monthly_9m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Monthly"] ~ 
                                    as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Monthly"]))
rslt_strat_monthly_9m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_monthly_9m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_strat_monthly_9m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Monthly" & kw_horizons$horizon ==1], probs=.75),
                               quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Monthly" & kw_horizons$horizon ==7], probs=.75),
                               quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Monthly" & kw_horizons$horizon ==35], probs=.75))

#kruskal wallis and dunn tests for 35 day ahead stratified rmse across depths and DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Daily"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Daily"]))
dunn_strat_daily_35day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Daily"] ~ 
                                     as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Daily"]))
rslt_strat_daily_35day=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_daily_35day$res, threshold = 0.05)$Letter)
rslt_strat_daily_35day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Daily" & kw_horizons$depth ==1], probs=.75),
                                quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Daily" & kw_horizons$depth ==5], probs=.75),
                                quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Daily" & kw_horizons$depth ==9], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Weekly"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Weekly"]))
dunn_strat_weekly_35day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Weekly"] ~ 
                                      as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Weekly"]))
rslt_strat_weekly_35day=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_weekly_35day$res, threshold = 0.05)$Letter)
rslt_strat_weekly_35day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Weekly" & kw_horizons$depth ==1], probs=.75),
                                 quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Weekly" & kw_horizons$depth ==5], probs=.75),
                                 quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Weekly" & kw_horizons$depth ==9], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Fortnightly"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Fortnightly"]))
dunn_strat_fortnightly_35day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Fortnightly"] ~ 
                                           as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Fortnightly"]))
rslt_strat_fortnightly_35day=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_fortnightly_35day$res, threshold = 0.05)$Letter)
rslt_strat_fortnightly_35day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Fortnightly" & kw_horizons$depth ==1], probs=.75),
                                      quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Fortnightly" & kw_horizons$depth ==5], probs=.75),
                                      quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Fortnightly" & kw_horizons$depth ==9], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Monthly"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Monthly"]))
dunn_strat_monthly_35day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Monthly"] ~ 
                                       as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Monthly"]))
rslt_strat_monthly_35day=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_monthly_35day$res, threshold = 0.05)$Letter)
rslt_strat_monthly_35day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Monthly" & kw_horizons$depth ==1], probs=.75),
                                  quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Monthly" & kw_horizons$depth ==5], probs=.75),
                                  quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35 & kw_horizons$model_id == "Monthly" & kw_horizons$depth ==9], probs=.75))



#kruskal wallis and dunn tests for 1m mixed rmse across horizon and DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Daily"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Daily"]))
dunn_mix_daily_1m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Daily"] ~ 
                                as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Daily"]))
rslt_mix_daily_1m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_daily_1m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_mix_daily_1m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==1], probs=.75),
                           quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==7], probs=.75),
                           quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly"]))
dunn_mix_weekly_1m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly"] ~ 
                                 as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly"]))
rslt_mix_weekly_1m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_weekly_1m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_mix_weekly_1m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==1], probs=.75),
                            quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==7], probs=.75),
                            quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly"]))
dunn_mix_fortnightly_1m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly"] ~ 
                                      as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly"]))
rslt_mix_fortnightly_1m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_fortnightly_1m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_mix_fortnightly_1m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==1], probs=.75),
                                 quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==7], probs=.75),
                                 quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly"]))
dunn_mix_monthly_1m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly"] ~ 
                                  as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly"]))
rslt_mix_monthly_1m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_monthly_1m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_mix_monthly_1m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly" & kw_horizons$horizon ==1], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly" & kw_horizons$horizon ==7], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly" & kw_horizons$horizon ==35], probs=.75))

#kruskal wallis and dunn tests for 1day mixed rmse across depths and DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Daily"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Daily"]))
dunn_mix_daily_1day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Daily"] ~ 
                                  as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Daily"]))
rslt_mix_daily_1day=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_daily_1day$res, threshold = 0.05)$Letter)
rslt_mix_daily_1day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Daily" & kw_horizons$depth ==1], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Daily" & kw_horizons$depth ==5], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Daily" & kw_horizons$depth ==9], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Weekly"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Weekly"]))
dunn_mix_weekly_1day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Weekly"] ~ 
                                   as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Weekly"]))
rslt_mix_weekly_1day=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_weekly_1day$res, threshold = 0.05)$Letter)
rslt_mix_weekly_1day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Weekly" & kw_horizons$depth ==1], probs=.75),
                              quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Weekly" & kw_horizons$depth ==5], probs=.75),
                              quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Weekly" & kw_horizons$depth ==9], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Fortnightly"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Fortnightly"]))
dunn_mix_fortnightly_1day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Fortnightly"] ~ 
                                        as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Fortnightly"]))
rslt_mix_fortnightly_1day=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_fortnightly_1day$res, threshold = 0.05)$Letter)
rslt_mix_fortnightly_1day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Fortnightly" & kw_horizons$depth ==1], probs=.75),
                                   quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Fortnightly" & kw_horizons$depth ==5], probs=.75),
                                   quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Fortnightly" & kw_horizons$depth ==9], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Monthly"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Monthly"]))
dunn_mix_monthly_1day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Monthly"] ~ 
                                    as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Monthly"]))
rslt_mix_monthly_1day=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_monthly_1day$res, threshold = 0.05)$Letter)
rslt_mix_monthly_1day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Monthly" & kw_horizons$depth ==1], probs=.75),
                               quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Monthly" & kw_horizons$depth ==5], probs=.75),
                               quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1 & kw_horizons$model_id == "Monthly" & kw_horizons$depth ==9], probs=.75))


#kruskal wallis and dunn tests for 5m mixed rmse across horizon and DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Daily"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Daily"]))
dunn_mix_daily_5m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Daily"] ~ 
                                as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Daily"]))
rslt_mix_daily_5m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_daily_5m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_mix_daily_5m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==1], probs=.75),
                           quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==7], probs=.75),
                           quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly"]))
dunn_mix_weekly_5m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly"] ~ 
                                 as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly"]))
rslt_mix_weekly_5m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_weekly_5m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_mix_weekly_5m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==1], probs=.75),
                            quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==7], probs=.75),
                            quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly"]))
dunn_mix_fortnightly_5m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly"] ~ 
                                      as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly"]))
rslt_mix_fortnightly_5m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_fortnightly_5m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_mix_fortnightly_5m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==1], probs=.75),
                                 quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==7], probs=.75),
                                 quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$v == "Fortnightly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Monthly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Monthly"]))
dunn_mix_monthly_5m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Monthly"] ~ 
                                  as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Monthly"]))
rslt_mix_monthly_5m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_monthly_5m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_mix_monthly_5m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Monthly" & kw_horizons$horizon ==1], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Monthly" & kw_horizons$horizon ==7], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Monthly" & kw_horizons$horizon ==35], probs=.75))

#kruskal wallis and dunn tests for 7day ahead mixed rmse across depths and DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Daily"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Daily"]))
dunn_mix_daily_7day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Daily"] ~ 
                                  as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Daily"]))
rslt_mix_daily_7day=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_daily_7day$res, threshold = 0.05)$Letter)
rslt_mix_daily_7day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Daily" & kw_horizons$depth ==1], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Daily" & kw_horizons$depth ==5], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Daily" & kw_horizons$depth ==9], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Weekly"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Weekly"]))
dunn_mix_weekly_7day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Weekly"] ~ 
                                   as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Weekly"]))
rslt_mix_weekly_7day=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_weekly_7day$res, threshold = 0.05)$Letter)
rslt_mix_weekly_7day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Weekly" & kw_horizons$depth ==1], probs=.75),
                              quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Weekly" & kw_horizons$depth ==5], probs=.75),
                              quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Weekly" & kw_horizons$depth ==9], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Fortnightly"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Fortnightly"]))
dunn_mix_fortnightly_7day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Fortnightly"] ~ 
                                        as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Fortnightly"]))
rslt_mix_fortnightly_7day=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_fortnightly_7day$res, threshold = 0.05)$Letter)
rslt_mix_fortnightly_7day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Fortnightly" & kw_horizons$depth ==1], probs=.75),
                                   quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Fortnightly" & kw_horizons$depth ==5], probs=.75),
                                   quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Fortnightly" & kw_horizons$depth ==9], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Monthly"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Monthly"]))
dunn_mix_monthly_7day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Monthly"] ~ 
                                    as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Monthly"]))
rslt_mix_monthly_7day=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_monthly_7day$res, threshold = 0.05)$Letter)
rslt_mix_monthly_7day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Monthly" & kw_horizons$depth ==1], probs=.75),
                               quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Monthly" & kw_horizons$depth ==5], probs=.75),
                               quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7 & kw_horizons$model_id == "Monthly" & kw_horizons$depth ==9], probs=.75))


#kruskal wallis and dunn tests for 9m mixed rmse across horizon and DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Daily"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Daily"]))
dunn_mix_daily_9m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Daily"] ~ 
                                as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Daily"]))
rslt_mix_daily_9m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_daily_9m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_mix_daily_9m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==1], probs=.75),
                           quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==7], probs=.75),
                           quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly"]))
dunn_mix_weekly_9m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly"] ~ 
                                 as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly"]))
rslt_mix_weekly_9m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_weekly_9m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_mix_weekly_9m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==1], probs=.75),
                            quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==7], probs=.75),
                            quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly"]))
dunn_mix_fortnightly_9m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly"] ~ 
                                      as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly"]))
rslt_mix_fortnightly_9m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_fortnightly_9m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_mix_fortnightly_9m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==1], probs=.75),
                                 quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==7], probs=.75),
                                 quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Monthly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Monthly"]))
dunn_mix_monthly_9m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Monthly"] ~ 
                                  as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Monthly"]))
rslt_mix_monthly_9m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_monthly_9m$res, threshold = 0.05)$Letter[c(1,3,2)])
rslt_mix_monthly_9m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Monthly" & kw_horizons$horizon ==1], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Monthly" & kw_horizons$horizon ==7], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Monthly" & kw_horizons$horizon ==35], probs=.75))

#kruskal wallis and dunn tests for 35 day ahead mixed rmse across depths and DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Daily"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Daily"]))
dunn_mix_daily_35day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Daily"] ~ 
                                   as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Daily"]))
rslt_mix_daily_35day=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_daily_35day$res, threshold = 0.05)$Letter)
rslt_mix_daily_35day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Daily" & kw_horizons$depth ==1], probs=.75),
                              quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Daily" & kw_horizons$depth ==5], probs=.75),
                              quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Daily" & kw_horizons$depth ==9], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Weekly"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Weekly"]))
dunn_mix_weekly_35day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Weekly"] ~ 
                                    as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Weekly"]))
rslt_mix_weekly_35day=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_weekly_35day$res, threshold = 0.05)$Letter)
rslt_mix_weekly_35day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Weekly" & kw_horizons$depth ==1], probs=.75),
                               quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Weekly" & kw_horizons$depth ==5], probs=.75),
                               quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Weekly" & kw_horizons$depth ==9], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Fortnightly"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Fortnightly"]))
dunn_mix_fortnightly_35day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Fortnightly"] ~ 
                                         as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Fortnightly"]))
rslt_mix_fortnightly_35day=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_fortnightly_35day$res, threshold = 0.05)$Letter)
rslt_mix_fortnightly_35day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Fortnightly" & kw_horizons$depth ==1], probs=.75),
                                    quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Fortnightly" & kw_horizons$depth ==5], probs=.75),
                                    quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Fortnightly" & kw_horizons$depth ==9], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Monthly"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Monthly"]))
dunn_mix_monthly_35day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Monthly"] ~ 
                                     as.factor(kw_horizons$depth[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Monthly"]))
rslt_mix_monthly_35day=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_monthly_35day$res, threshold = 0.05)$Letter)
rslt_mix_monthly_35day_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Monthly" & kw_horizons$depth ==1], probs=.75),
                                quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Monthly" & kw_horizons$depth ==5], probs=.75),
                                quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35 & kw_horizons$model_id == "Monthly" & kw_horizons$depth ==9], probs=.75))


#create table to export mixed and stratified p-vals comparing horizons for each depth and DA freq
pvals <- data.frame("Depth_m"= c(rep(1,12),rep(5,12),rep(9,12)), 
                    "DA" = rep(c(rep("daily",3),rep("weekly",3),rep("fortnightly",3),rep("monthly",3)),3),
                    "Comparison" = rep(c(dunn_mix_daily_1m$res$Comparison),12),
                    # "Z" = c(dunn_mix_daily_1m$res$Z,dunn_mix_weekly_1m$res$Z,dunn_mix_fortnightly_1m$res$Z,dunn_mix_monthly_1m$res$Z,
                    #         dunn_mix_daily_5m$res$Z,dunn_mix_weekly_5m$res$Z,dunn_mix_fortnightly_5m$res$Z,dunn_mix_monthly_5m$res$Z,
                    #         dunn_mix_daily_9m$res$Z,dunn_mix_weekly_9m$res$Z,dunn_mix_fortnightly_9m$res$Z,dunn_mix_monthly_9m$res$Z),
                    "Mixed_pvalue" = c(dunn_mix_daily_1m$res$P.adj,dunn_mix_weekly_1m$res$P.adj,dunn_mix_fortnightly_1m$res$P.adj,dunn_mix_monthly_1m$res$P.adj,
                                       dunn_mix_daily_5m$res$P.adj,dunn_mix_weekly_5m$res$P.adj,dunn_mix_fortnightly_5m$res$P.adj,dunn_mix_monthly_5m$res$P.adj,
                                       dunn_mix_daily_9m$res$P.adj,dunn_mix_weekly_9m$res$P.adj,dunn_mix_fortnightly_9m$res$P.adj,dunn_mix_monthly_9m$res$P.adj),
                    "Stratified_pvalue" = c(dunn_strat_daily_1m$res$P.adj,dunn_strat_weekly_1m$res$P.adj,dunn_strat_fortnightly_1m$res$P.adj,dunn_strat_monthly_1m$res$P.adj,
                                            dunn_strat_daily_5m$res$P.adj,dunn_strat_weekly_5m$res$P.adj,dunn_strat_fortnightly_5m$res$P.adj,dunn_strat_monthly_5m$res$P.adj,
                                            dunn_strat_daily_9m$res$P.adj,dunn_strat_weekly_9m$res$P.adj,dunn_strat_fortnightly_9m$res$P.adj,dunn_strat_monthly_9m$res$P.adj))
#write.csv(pvals,file.path(lake_directory,"analysis/data/UC_pvals_horizon.csv"),row.names = FALSE)

#create table to export mixed and stratified p-vals comparing depths for each horizon and DA freq
pvals_depths <- data.frame("Horizon_days"= c(rep(1,12),rep(7,12),rep(35,12)), 
                           "DA" = rep(c(rep("daily",3),rep("weekly",3),rep("fortnightly",3),rep("monthly",3)),3),
                           "Comparison" = rep(c(dunn_mix_daily_1day$res$Comparison),12),
                           "Mixed_pvalue" = c(dunn_mix_daily_1day$res$P.adj,dunn_mix_weekly_1day$res$P.adj,dunn_mix_fortnightly_1day$res$P.adj,dunn_mix_monthly_1day$res$P.adj,
                                              dunn_mix_daily_7day$res$P.adj,dunn_mix_weekly_7day$res$P.adj,dunn_mix_fortnightly_7day$res$P.adj,dunn_mix_monthly_7day$res$P.adj,
                                              dunn_mix_daily_35day$res$P.adj,dunn_mix_weekly_35day$res$P.adj,dunn_mix_fortnightly_35day$res$P.adj,dunn_mix_monthly_35day$res$P.adj),
                           "Stratified_pvalue" = c(dunn_strat_daily_1day$res$P.adj,dunn_strat_weekly_1day$res$P.adj,dunn_strat_fortnightly_1day$res$P.adj,dunn_strat_monthly_1day$res$P.adj,
                                                   dunn_strat_daily_7day$res$P.adj,dunn_strat_weekly_7day$res$P.adj,dunn_strat_fortnightly_7day$res$P.adj,dunn_strat_monthly_7day$res$P.adj,
                                                   dunn_strat_daily_35day$res$P.adj,dunn_strat_weekly_35day$res$P.adj,dunn_strat_fortnightly_35day$res$P.adj,dunn_strat_monthly_35day$res$P.adj))
#write.csv(pvals_depths,file.path(lake_directory,"analysis/data/UC_pvals_depth.csv"),row.names = FALSE)


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#### FIGURE 5: 1,5,9m depth facets for each DA frequency for 1,7,and 35-day ahead forecasts  ####

letters <- data.frame("depth"= c(rep(1,8),rep(5,8),rep(9,8)),
                      "horizon" = rep(c(1,7,35),8),
                      "model_id" = rep(c("Daily","Weekly","Fortnightly","Monthly"),6),
                      "x" = rep(c(1,2,3,4),6),
                      "phen" = rep(c(rep("Mixed",4),rep("Stratified",4)),3),
                      "letter" = tolower(c(rslt_mix_1m, rslt_strat_1m, rslt_mix_5m, rslt_strat_5m,
                                           rslt_mix_9m,rslt_strat_9m)),
                      "max.RMSE" = c(rep(2,16),rep(2,8)))

#rename depth facets
depths <- c("1m","5m","9m")
names(depths) <- c("1","5","9")

#order factor levels
kw_horizons$model_id <- factor(kw_horizons$model_id, levels = c("Daily", "Weekly", "Fortnightly", "Monthly"))

ggplot(kw_horizons, aes(model_id, RMSE, fill=as.factor(horizon))) +  ylab("RMSE") + xlab("")+
  geom_bar(stat="identity",position="dodge") + theme_bw() + guides(fill=guide_legend(title="Horizon")) +
  geom_hline(yintercept=2, linetype='dashed', col = 'black', size=0.3)+ 
  scale_fill_manual(values=c("#81A665","#E0CB48","#D08151"),labels = c("1day", "7day", "35day")) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = c(0.75,0.31),
        legend.background = element_blank(),legend.direction = "horizontal", panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0.05,-0.2,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 6),legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  geom_text(data=letters,aes(x=model_id,y=0.2+max.RMSE,label=letter),hjust=0.1,vjust = -0.1, size=2.5) +
  facet_grid(depth~phen, scales="free_y",labeller = labeller(depth = depths)) 
ggsave(file.path(lake_directory,"analysis/figures/UC_RMSEvsDAfreq_depth_facets_fig5.jpg"),width=3.5, height=4)


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#### FIGURE 6: mixed and stratified RMSE tileplots aggregated across all depths  ####

#round rmse to nearest 0.5 for tile plot below
forecast_horizon_avg$RMSE_bins <- plyr::round_any(forecast_horizon_avg$RMSE,0.5) 

#figure for horizon vs frequency to compare forecast skill
ggplot(forecast_horizon_avg, aes(model_id, horizon, fill=RMSE_bins)) + 
  ylab("Horizon (days)") + theme_bw() + facet_wrap(~phen) + geom_tile(width=0.8) + xlab("") +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.spacing=unit(0, "cm"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  guides(fill=guide_legend(title="RMSE")) +  scale_fill_gradientn(colors = hcl.colors(4, "BuPu")) 
ggsave(file.path(lake_directory,"analysis/figures/UC_HorizonvsDA_tileplot_fig6.jpg"),width=4, height=3.5)

#median mixed vs stratified period rmse across different horizons
median(forecast_horizon_avg$RMSE_bins[forecast_horizon_avg$phen=="Mixed" & forecast_horizon_avg$horizon==9 & forecast_horizon_avg$model_id=="Fortnightly"])
median(forecast_horizon_avg$RMSE_bins[forecast_horizon_avg$phen=="Stratified" & forecast_horizon_avg$horizon==1 & forecast_horizon_avg$model_id=="Fortnightly"])

#------------------------------------------------------------------------------------------------#
# Data assimilation figure 4
DA <- all_DA_forecasts

#pull out 2 horizons for fig 4 (mixed vs stratified)
DA_sub <- DA[(DA$datetime >="2021-01-01" & DA$datetime <="2021-01-31") | (DA$datetime >="2021-06-24" & DA$datetime <="2021-07-25"),]

#summary df to average the forecasts for each DA freq, horizon, depth, and date
DA_sub_final <- plyr::ddply(DA_sub, c("depth", "datetime","horizon", "model_id"), function(x) {
  data.frame(
    forecast_mean = mean(x$mean, na.rm = TRUE),
    forecast_sd = mean(x$sd, na.rm = TRUE),
    forecast_upper_95 = mean(x$quantile97.5, na.rm=TRUE),
    forecast_lower_95 = mean(x$quantile02.5, na.rm=TRUE),
    observed = mean(x$observation, na.rm=TRUE)
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

#change datetime format
DA_sub_final$datetime <- as.Date(DA_sub_final$datetime)

#add in phen column
DA_sub_final$phen <- "Stratified" 
DA_sub_final$phen[DA_sub_final$datetime >="2021-01-01" & DA_sub_final$date <="2021-01-31"] <- "Mixed"

#change DA factor order
DA_sub_final$model_id <- factor(DA_sub_final$model_id, levels = c("Daily", "Weekly","Fortnightly","Monthly"))

#change order of DA frequencies so daily is plotted on top
DA_sub_final$model_id <- factor(DA_sub_final$model_id, levels=rev(levels(DA_sub_final$model_id)))

ggplot(subset(DA_sub_final, horizon==1 & depth %in% c(1,5,9)), aes(datetime, forecast_mean, color=model_id)) + 
  geom_ribbon(aes(x=datetime, y = forecast_mean, ymin = forecast_mean-forecast_sd, ymax = forecast_mean+forecast_sd, 
                  color=model_id, fill=model_id),alpha=0.4)  + theme_bw() + 
  #geom_line(aes(datetime, forecast_mean), size=0.2, data = . %>% filter(model_id %in% c("Daily"))) +
  geom_point(aes(x=datetime, y=observed), col="black", size=0.3) +  
  facet_wrap(depth~phen, scales = "free", labeller = labeller(depth=depths,.multi_line = FALSE),ncol = 2) + scale_x_date() +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = c(0.76,0.985),
        legend.background = element_blank(),legend.direction = "horizontal", panel.grid.minor = element_blank(), legend.key=element_rect(fill=NA),
        plot.margin = unit(c(0,0.05,-0.2,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 4.5),legend.text  = element_text(size = 4.5), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  ylab(expression("Temperature ("*~degree*C*")")) + xlab("")  + scale_color_manual(values=rev(cb_friendly_2)) + scale_fill_manual(values=rev(cb_friendly_2)) +
  guides(fill = guide_legend(title="", override.aes = list(alpha=1),reverse = TRUE), color="none")
ggsave(file.path(lake_directory,"analysis/figures/UC_AssimilationVSdatafreq_fig4.jpg"), width=3.5, height=4) 

#uncertainty range across depths and mixed vs stratified periods
((mean(DA_sub_final$forecast_upper_95[DA_sub_final$phen=="Mixed" & DA_sub_final$model_id=="Monthly"]) -
    mean(DA_sub_final$forecast_lower_95[DA_sub_final$phen=="Mixed" & DA_sub_final$model_id=="Monthly"])) -
    
    (mean(DA_sub_final$forecast_upper_95[DA_sub_final$phen=="Mixed" & DA_sub_final$model_id=="Daily"]) -
       mean(DA_sub_final$forecast_lower_95[DA_sub_final$phen=="Mixed" & DA_sub_final$model_id=="Daily"]))) /
  
  mean((mean(DA_sub_final$forecast_upper_95[DA_sub_final$phen=="Mixed" & DA_sub_final$model_id=="Monthly"]) -
          mean(DA_sub_final$forecast_lower_95[DA_sub_final$phen=="Mixed" & DA_sub_final$model_id=="Monthly"])),
       (mean(DA_sub_final$forecast_upper_95[DA_sub_final$phen=="Mixed" & DA_sub_final$model_id=="Daily"]) -
          mean(DA_sub_final$forecast_lower_95[DA_sub_final$phen=="Mixed" & DA_sub_final$model_id=="Daily"]))) *100


((mean(DA_sub_final$forecast_upper_95[DA_sub_final$phen=="Stratified" & DA_sub_final$model_id=="Monthly"]) -
    mean(DA_sub_final$forecast_lower_95[DA_sub_final$phen=="Stratified" & DA_sub_final$model_id=="Monthly"])) -
    
    (mean(DA_sub_final$forecast_upper_95[DA_sub_final$phen=="Stratified" & DA_sub_final$model_id=="Daily"]) -
       mean(DA_sub_final$forecast_lower_95[DA_sub_final$phen=="Stratified" & DA_sub_final$model_id=="Daily"]))) / 
  
  mean((mean(DA_sub_final$forecast_upper_95[DA_sub_final$phen=="Stratified" & DA_sub_final$model_id=="Monthly"]) -
          mean(DA_sub_final$forecast_lower_95[DA_sub_final$phen=="Stratified" & DA_sub_final$model_id=="Monthly"])),
       (mean(DA_sub_final$forecast_upper_95[DA_sub_final$phen=="Stratified" & DA_sub_final$model_id=="Daily"]) -
          mean(DA_sub_final$forecast_lower_95[DA_sub_final$phen=="Stratified" & DA_sub_final$model_id=="Daily"]))) *100

#-------------------------------------------------------------------------------#
# Fig 7 - fig to compare forecast skill across horizons and depths 

#more KW tests
#kruskal wallis and dunn tests for 1m + 1day Mixed rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon==1] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon==1])
dunn_mix_1m_1d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon==1] ~ 
                             kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon==1])
rslt_mix_1m_1d=cldList(P.adj ~ Comparison, data=dunn_mix_1m_1d$res, threshold = 0.05)$Letter[c(1,4,2,3)]
median_mix_1m_1d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Daily"])
median_mix_1m_1d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Weekly"])
median_mix_1m_1d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Fortnightly"])
median_mix_1m_1d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 1m + 7day Mixed rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon==7] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon==7])
dunn_mix_1m_7d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon==7] ~ 
                             kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon==7])
rslt_mix_1m_7d=cldList(P.adj ~ Comparison, data=dunn_mix_1m_7d$res, threshold = 0.05)$Letter[c(1,4,2,3)]
median_mix_1m_7d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Daily"])
median_mix_1m_7d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Weekly"])
median_mix_1m_7d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Fortnightly"])
median_mix_1m_7d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 1m + 35day Mixed rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon==35] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon==35])
dunn_mix_1m_35d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon==35] ~ 
                              kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon==35])
rslt_mix_1m_35d=cldList(P.adj ~ Comparison, data=dunn_mix_1m_35d$res, threshold = 0.05)$Letter[c(1,4,2,3)]
median_mix_1m_35d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Daily"])
median_mix_1m_35d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Weekly"])
median_mix_1m_35d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Fortnightly"])
median_mix_1m_35d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 5m + 1day Mixed rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon==1] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon==1])
dunn_mix_5m_1d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon==1] ~ 
                             kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon==1])
rslt_mix_5m_1d=cldList(P.adj ~ Comparison, data=dunn_mix_5m_1d$res, threshold = 0.05)$Letter[c(1,4,2,3)]
median_mix_5m_1d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Daily"])
median_mix_5m_1d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Weekly"])
median_mix_5m_1d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Fortnightly"])
median_mix_5m_1d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 5m + 7day Mixed rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon==7] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon==7])
dunn_mix_5m_7d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon==7] ~ 
                             kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon==7])
rslt_mix_5m_7d=cldList(P.adj ~ Comparison, data=dunn_mix_5m_7d$res, threshold = 0.05)$Letter[c(1,4,2,3)]
median_mix_5m_7d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Daily"])
median_mix_5m_7d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Weekly"])
median_mix_5m_7d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Fortnightly"])
median_mix_5m_7d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 5m + 35day Mixed rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon==35] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon==35])
dunn_mix_5m_35d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon==35] ~ 
                              kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon==35])
rslt_mix_5m_35d=cldList(P.adj ~ Comparison, data=dunn_mix_5m_35d$res, threshold = 0.05)$Letter[c(1,4,2,3)]
median_mix_5m_35d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Daily"])
median_mix_5m_35d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Weekly"])
median_mix_5m_35d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Fortnightly"])
median_mix_5m_35d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 9m + 1day Mixed rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon==1] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon==1])
dunn_mix_9m_1d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon==1] ~ 
                             kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon==1])
rslt_mix_9m_1d=cldList(P.adj ~ Comparison, data=dunn_mix_9m_1d$res, threshold = 0.05)$Letter[c(1,4,2,3)]
median_mix_9m_1d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Daily"])
median_mix_9m_1d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Weekly"])
median_mix_9m_1d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Fortnightly"])
median_mix_9m_1d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 9m + 7day Mixed rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon==7] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon==7])
dunn_mix_9m_7d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon==7] ~ 
                             kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon==7])
rslt_mix_9m_7d=cldList(P.adj ~ Comparison, data=dunn_mix_9m_7d$res, threshold = 0.05)$Letter[c(1,4,2,3)]
median_mix_9m_7d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Daily"])
median_mix_9m_7d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Weekly"])
median_mix_9m_7d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Fortnightly"])
median_mix_9m_7d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 9m + 35day Mixed rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon==35] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon==35])
dunn_mix_9m_35d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon==35] ~ 
                              kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon==35])
rslt_mix_9m_35d=cldList(P.adj ~ Comparison, data=dunn_mix_9m_35d$res, threshold = 0.05)$Letter[c(1,4,2,3)]
median_mix_9m_35d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Daily"])
median_mix_9m_35d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Weekly"])
median_mix_9m_35d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Fortnightly"])
median_mix_9m_35d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 1m + 1day stratified rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon==1] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon==1])
dunn_strat_1m_1d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon==1] ~ 
                               kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon==1])
rslt_strat_1m_1d=cldList(P.adj ~ Comparison, data=dunn_strat_1m_1d$res, threshold = 0.05)$Letter[c(1,4,2,3)]
median_strat_1m_1d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Daily"])
median_strat_1m_1d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Weekly"])
median_strat_1m_1d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Fortnightly"])
median_strat_1m_1d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 1m + 7day stratified rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon==7] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon==7])
dunn_strat_1m_7d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon==7] ~ 
                               kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon==7])
rslt_strat_1m_7d=cldList(P.adj ~ Comparison, data=dunn_strat_1m_7d$res, threshold = 0.05)$Letter[c(1,4,2,3)]
median_strat_1m_7d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Daily"])
median_strat_1m_7d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Weekly"])
median_strat_1m_7d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Fortnightly"])
median_strat_1m_7d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 1m + 35day stratified rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon==35] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon==35])
dunn_strat_1m_35d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon==35] ~ 
                                kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon==35])
rslt_strat_1m_35d=cldList(P.adj ~ Comparison, data=dunn_strat_1m_35d$res, threshold = 0.05)$Letter[c(1,4,2,3)]
median_strat_1m_35d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Daily"])
median_strat_1m_35d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Weekly"])
median_strat_1m_35d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Fortnightly"])
median_strat_1m_35d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 5m + 1day stratified rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon==1] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon==1])
dunn_strat_5m_1d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon==1] ~ 
                               kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon==1])
rslt_strat_5m_1d=cldList(P.adj ~ Comparison, data=dunn_strat_5m_1d$res, threshold = 0.05)$Letter[c(1,4,2,3)]
median_strat_5m_1d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Daily"])
median_strat_5m_1d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Weekly"])
median_strat_5m_1d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Fortnightly"])
median_strat_5m_1d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 5m + 7day stratified rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon==7] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon==7])
dunn_strat_5m_7d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon==7] ~ 
                               kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon==7])
rslt_strat_5m_7d=cldList(P.adj ~ Comparison, data=dunn_strat_5m_7d$res, threshold = 0.05)$Letter[c(1,4,2,3)]
median_strat_5m_7d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Daily"])
median_strat_5m_7d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Weekly"])
median_strat_5m_7d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Fortnightly"])
median_strat_5m_7d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 5m + 35day stratified rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon==35] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon==35])
dunn_strat_5m_35d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon==35] ~ 
                                kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon==35])
rslt_strat_5m_35d=cldList(P.adj ~ Comparison, data=dunn_strat_5m_35d$res, threshold = 0.05)$Letter[c(1,4,2,3)]
median_strat_5m_35d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Daily"])
median_strat_5m_35d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Weekly"])
median_strat_5m_35d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Fortnightly"])
median_strat_5m_35d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 9m + 1day stratified rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon==1] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon==1])
dunn_strat_9m_1d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon==1] ~ 
                               kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon==1])
rslt_strat_9m_1d=cldList(P.adj ~ Comparison, data=dunn_strat_9m_1d$res, threshold = 0.05)$Letter[c(1,4,2,3)]
median_strat_9m_1d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Daily"])
median_strat_9m_1d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Weekly"])
median_strat_9m_1d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Fortnightly"])
median_strat_9m_1d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==1 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 9m + 7day stratified rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon==7] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon==7])
dunn_strat_9m_7d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon==7] ~ 
                               kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon==7])
rslt_strat_9m_7d=cldList(P.adj ~ Comparison, data=dunn_strat_9m_7d$res, threshold = 0.05)$Letter[c(1,4,2,3)]
median_strat_9m_7d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Daily"])
median_strat_9m_7d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Weekly"])
median_strat_9m_7d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Fortnightly"])
median_strat_9m_7d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==7 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 9m + 35day stratified rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon==35] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon==35])
dunn_strat_9m_35d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon==35] ~ 
                                kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon==35])
rslt_strat_9m_35d=cldList(P.adj ~ Comparison, data=dunn_strat_9m_35d$res, threshold = 0.05)$Letter[c(1,4,2,3)]
median_strat_9m_35d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Daily"])
median_strat_9m_35d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Weekly"])
median_strat_9m_35d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Fortnightly"])
median_strat_9m_35d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$horizon ==35 & kw_horizons$model_id=="Monthly"])

#create table to export mixed and stratified p-vals comparing depths for each horizon and DA freq
pvals_fig7_DA <- data.frame("Depth_m"= c(rep(1,18),rep(5,18),rep(9,18)),
                            "Horizon" = c(rep(c(rep(1,6),rep(7,6),rep(35,6)),3)),
                            "Comparison" = rep(c(dunn_mix_1m_1d$res$Comparison),9),
                            "Mixed_pvalue" = c(dunn_mix_1m_1d$res$P.adj,dunn_mix_1m_7d$res$P.adj,dunn_mix_1m_35d$res$P.adj,
                                               dunn_mix_5m_1d$res$P.adj,dunn_mix_5m_7d$res$P.adj,dunn_mix_5m_35d$res$P.adj,
                                               dunn_mix_9m_1d$res$P.adj,dunn_mix_9m_7d$res$P.adj,dunn_mix_9m_35d$res$P.adj),
                            "Stratified_pvalue" = c(dunn_strat_1m_1d$res$P.adj,dunn_strat_1m_7d$res$P.adj,dunn_strat_1m_35d$res$P.adj,
                                                    dunn_strat_5m_1d$res$P.adj,dunn_strat_5m_7d$res$P.adj,dunn_strat_5m_35d$res$P.adj,
                                                    dunn_strat_9m_1d$res$P.adj,dunn_strat_9m_7d$res$P.adj,dunn_strat_9m_35d$res$P.adj))
#write.csv(pvals_fig7_DA,file.path(lake_directory,"analysis/data/UC_pvals_fig7_DA.csv"),row.names = FALSE)



#create new df with all letters 
DA_horizon_letters <- data.frame("Depth_m"= c(rep(c(rep(1,12),rep(5,12),rep(9,12)),2)),
                                 "Horizon" = c(rep(c(1,7,35),6)),
                                 "TempDynamics" = c(rep("Mixed",36),rep("Stratified",36)),
                                 "model_id" = c(rep(c("Daily","Weekly","Fortnightly","Monthly"),18)),
                                 "letter" = tolower(c(rslt_mix_1m_1d,rslt_mix_1m_7d,rslt_mix_1m_35d,
                                                      rslt_mix_5m_1d,rslt_mix_5m_7d,rslt_mix_5m_35d,
                                                      rslt_mix_9m_1d,rslt_mix_9m_7d,rslt_mix_9m_35d,
                                                      rslt_strat_1m_1d,rslt_strat_1m_7d,rslt_strat_1m_35d,
                                                      rslt_strat_5m_1d,rslt_strat_5m_7d,rslt_strat_5m_35d,
                                                      rslt_strat_9m_1d,rslt_strat_9m_7d,rslt_strat_9m_35d)))

#change factor order
DA_horizon_letters$model_id <- factor(DA_horizon_letters$model_id, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))

#now make a figure...
ggplot(DA_horizon_letters, aes(model_id, as.factor(Horizon), fill=as.factor(letter))) + 
  ylab("Horizon (days)") + xlab("") + theme_bw() + geom_tile(width=0.8) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = "right",
        legend.background = element_blank(),legend.direction = "vertical", panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0.05,-0.2,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 6),legend.text  = element_text(size = 8), panel.spacing=unit(0, "cm"), legend.margin=margin(0,0,0,0),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  facet_grid(Depth_m~TempDynamics, scales="free_y",labeller = labeller(Depth_m = depths)) +
  guides(fill=guide_legend(title="")) + scale_fill_brewer(type="div", palette=6, direction = 1)
#ggsave(file.path(lake_directory,"analysis/figures/UC_posthoc_letters_horizon_depth_facets_fig7.jpg"),width=3.5, height=4)

#------------------------------------------------------------------------------------------------#
#parameter evolution figs 
#read in all forecasts 
params_dir <- arrow::SubTreeFileSystem$create(file.path(lake_directory,"scores/UC"))
params <- arrow::open_dataset(score_dir) |> collect() |>   
  filter(variable %in% c("lw_factor","zone1temp","zone2temp"), horizon >=0)

#need to round horizon becuase they are in decimal form for some reason...
params$horizon <- ceiling(params$horizon)

#change datetime format
params$datetime <- as.Date(params$datetime)

#change model_id to be all uppercase
params$model_id <- str_to_title(params$model_id)

#change DA factor order
params$model_id <- factor(params$model_id, levels = c("Daily", "Weekly","Fortnightly","Monthly"))

ggplot(subset(params,horizon %in% c(1)), aes(datetime, mean, color=model_id)) + theme_bw() +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = c(0.6,0.98),
        legend.background = element_blank(),legend.direction = "horizontal", panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0.05,-0.2,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 6),legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  scale_color_manual(values=cb_friendly_2) +
  scale_fill_manual(values=cb_friendly_2) +
  facet_wrap(~variable, scales="free_y",ncol=1) + ylab("parameter")+ xlab("")+
  scale_x_date(date_labels = "%b") + #ylab(expression("Temperature ("*~degree*C*")")) 
  geom_ribbon(aes(y = mean, ymin = mean-sd, ymax = mean+sd, color=model_id, fill=model_id), alpha=0.5) +
  guides(fill = guide_legend(title="DA frequency", override.aes = list(alpha=1)), color="none")
#ggsave(file.path(lake_directory,"analysis/figures/UC_paramevolvsHorizon.jpg"),width=3.5, height=4)

#figuring out the date that DA parameters diverge
mean(params$mean[params$variable=="lw_factor" & params$model_id=="Daily" & params$datetime >= "2021-04-01"])

mean(params$mean[params$variable=="lw_factor" & params$model_id=="Weekly" & params$datetime >= "2021-04-01"])
mean(params$mean[params$variable=="lw_factor" & params$model_id=="Fortnightly" & params$datetime >= "2021-04-01"])
mean(params$mean[params$variable=="lw_factor" & params$model_id=="Monthly" & params$datetime >= "2021-04-01"])

mean(c(last(params$mean[params$variable=="lw_factor" & params$model_id=="Weekly"]),
       last(params$mean[params$variable=="lw_factor" & params$model_id=="Fortnightly"]),
       last(params$mean[params$variable=="lw_factor" & params$model_id=="Monthly"])))

mean(c(last(params$mean[params$variable=="zone1temp" & params$model_id=="Daily"]),
       last(params$mean[params$variable=="zone1temp" & params$model_id=="Weekly"]),
       last(params$mean[params$variable=="zone1temp" & params$model_id=="Fortnightly"]),
       last(params$mean[params$variable=="zone1temp" & params$model_id=="Monthly"])))

mean(c(last(params$mean[params$variable=="zone2temp" & params$model_id=="Daily"]),
       last(params$mean[params$variable=="zone2temp" & params$model_id=="Weekly"]),
       last(params$mean[params$variable=="zone2temp" & params$model_id=="Fortnightly"]),
       last(params$mean[params$variable=="zone2temp" & params$model_id=="Monthly"])))

#------------------------------------------------------------------------------------------------#
# Forecasts with IC on vs off

#read in all forecasts with IC on
score_dir_yesIC <- arrow::SubTreeFileSystem$create(file.path(lake_directory,"scores/DA_study"))
all_DA_forecasts_yesIC <- arrow::open_dataset(score_dir_yesIC) |> collect() |>   
  filter(!is.na(observation), variable == "temperature",horizon > 0)

#round depths up to nearest m 
all_DA_forecasts_yesIC$depth <- ceiling(all_DA_forecasts_yesIC$depth)

#add a group number so that I can average horizons later on
all_DA_forecasts_yesIC <- all_DA_forecasts_yesIC %>% 
  mutate(group = case_when(all_DA_forecasts_yesIC$horizon <= 5 ~ "1-5",
                           all_DA_forecasts_yesIC$horizon <=10 & all_DA_forecasts_yesIC$horizon > 5 ~ "6-10",
                           all_DA_forecasts_yesIC$horizon <=15 & all_DA_forecasts_yesIC$horizon > 10 ~ "11-15",
                           all_DA_forecasts_yesIC$horizon <=20 & all_DA_forecasts_yesIC$horizon > 15 ~ "16-20",
                           all_DA_forecasts_yesIC$horizon <=25 & all_DA_forecasts_yesIC$horizon > 20 ~ "21-25",
                           all_DA_forecasts_yesIC$horizon <=30 & all_DA_forecasts_yesIC$horizon > 25 ~ "26-30",
                           all_DA_forecasts_yesIC$horizon <=36 & all_DA_forecasts_yesIC$horizon > 30 ~ "31-35"))

strat_date<- "2021-11-07"

#add stratified vs mixed col
all_DA_forecasts_yesIC$phen <- ifelse(all_DA_forecasts_yesIC$datetime <= as.POSIXct(strat_date) & 
                                        all_DA_forecasts_yesIC$datetime >="2021-03-13","Stratified", "Mixed")

#remove n=6 days with ice-cover 
all_DA_forecasts_yesIC <- all_DA_forecasts_yesIC[!(as.Date(all_DA_forecasts_yesIC$datetime) %in% c(as.Date("2021-01-10"), as.Date("2021-01-11"),as.Date("2021-01-30"),
                                                                                                   as.Date("2021-02-13"),as.Date("2021-02-14"),as.Date("2021-02-15"))),]

#change model_id to be all uppercase
all_DA_forecasts_yesIC$model_id <- str_to_title(all_DA_forecasts_yesIC$model_id)

#only keep 2021 data
all_DA_forecasts_yesIC <- all_DA_forecasts_yesIC[all_DA_forecasts_yesIC$datetime<="2021-12-31",]

#------------------------------------------------------------------------------#
#calculate forecast skill metrics

#forecast skill for each depth and horizon
forecast_skill_depth_horizon_yesIC <-  plyr::ddply(all_DA_forecasts_yesIC, c("depth","horizon","phen", "model_id"), function(x) {
  data.frame(
    RMSE = sqrt(mean((x$mean - x$observation)^2, na.rm = TRUE)),
    MAE = mean(abs(x$mean - x$observation), na.rm = TRUE),
    pbias = 100 * (sum(x$mean - x$observation, na.rm = TRUE) / sum(x$observation, na.rm = TRUE)),
    CRPS = verification::crps(x$observation, as.matrix(x[, c(7,9)]))$CRPS,
    variance = (mean(x$sd))^2
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 


#order DA frequencies
forecast_skill_depth_horizon_yesIC$model_id <- factor(forecast_skill_depth_horizon_yesIC$model_id, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))

#df with averaged forecast skill for all days (group by horizon, DA, and phen)
forecast_horizon_avg_yesIC <- plyr::ddply(all_DA_forecasts_yesIC, c("horizon", "model_id", "phen"), function(x) {
  data.frame(
    RMSE = sqrt(mean((x$mean - x$observation)^2, na.rm = TRUE)),
    MAE = mean(abs(x$mean - x$observation), na.rm = TRUE),
    pbias = 100 * (sum(x$mean - x$observation, na.rm = TRUE) / sum(x$observation, na.rm = TRUE)),
    CRPS = verification::crps(x$observation, as.matrix(x[, c(7,9)]))$CRPS
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

#order DA frequencies
forecast_horizon_avg_yesIC$model_id <- factor(forecast_horizon_avg_yesIC$model_id, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))

#df with averaged forecast skill for all days (group by depth, DA, and phen)
forecast_depth_avg_yesIC <- plyr::ddply(all_DA_forecasts_yesIC, c("depth", "model_id", "phen"), function(x) {
  data.frame(
    RMSE = sqrt(mean((x$mean - x$observation)^2, na.rm = TRUE)),
    MAE = mean(abs(x$mean - x$observation), na.rm = TRUE),
    pbias = 100 * (sum(x$mean - x$observation, na.rm = TRUE) / sum(x$observation, na.rm = TRUE)),
    CRPS = verification::crps(x$observation, as.matrix(x[, c(7,9)]))$CRPS,
    variance = (mean(x$sd))^2
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

#order DA frequencies
forecast_depth_avg_yesIC$model_id <- factor(forecast_depth_avg_yesIC$model_id, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))

#add IC y/n column
forecast_skill_depth_horizon$IC <- "no"
forecast_skill_depth_horizon_yesIC$IC <- "yes"

forecast_horizon_avg$IC <- "n"
forecast_horizon_avg_yesIC$IC <- "y"

forecast_depth_avg$IC <- "no"
forecast_depth_avg_yesIC$IC <- "yes"

#combine to make a massive df
UC <- rbind(forecast_skill_depth_horizon_yesIC,forecast_skill_depth_horizon)
UC_depth <- rbind(forecast_depth_avg_yesIC,forecast_depth_avg)

#round depth to nearest integer
UC$depth <- ceiling(UC$depth)
UC_depth$depth <- ceiling(UC_depth$depth)

ggplot(subset(UC_depth, depth %in% c(1,5,9)), aes(model_id, RMSE, fill=as.factor(IC))) +  ylab("RMSE") + xlab("")+
  geom_bar(stat="identity",position="dodge") + theme_bw() + guides(fill=guide_legend(title="")) + geom_hline(yintercept=2,linetype="dashed") +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = c(0.75,0.31),
        legend.background = element_blank(),legend.direction = "horizontal", panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0.05,-0.2,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 6),legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  facet_grid(depth~phen, scales="free",labeller = labeller(depth = depths)) + scale_fill_manual(values=c("#81A665","#E0CB48")) 
ggsave(file.path(lake_directory,"analysis/figures/UC_RMSEvsDAfreq_depth_facets_IC_allhorizons.jpg"),width=3.5, height=4)

ggplot(subset(UC, depth %in% c(1,5,9) & horizon==1), aes(model_id, RMSE, fill=as.factor(IC))) +  ylab("RMSE") + xlab("")+
  geom_bar(stat="identity",position="dodge") + theme_bw() + guides(fill=guide_legend(title="")) + geom_hline(yintercept=2,linetype="dashed") +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = c(0.75,0.31),
        legend.background = element_blank(),legend.direction = "horizontal", panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0.05,-0.2,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 6),legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  facet_grid(depth~phen, scales="free",labeller = labeller(depth = depths)) + scale_fill_manual(values=c("#81A665","#E0CB48")) 
ggsave(file.path(lake_directory,"analysis/figures/UC_RMSEvsDAfreq_depth_facets_IC_1day.jpg"),width=3.5, height=4)


ggplot(subset(UC_depth, depth %in% c(1,5,9)) ,aes(model_id, variance, fill=as.factor(IC))) +  ylab("variance") + xlab("")+
  geom_bar(stat="identity",position="dodge") + theme_bw() + guides(fill=guide_legend(title="")) + 
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = c(0.75,0.31),
        legend.background = element_blank(),legend.direction = "horizontal", panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0.05,-0.2,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 6),legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  facet_grid(depth~phen, scales="free",labeller = labeller(depth = depths)) + scale_fill_manual(values=c("#81A665","#E0CB48")) 
ggsave(file.path(lake_directory,"analysis/figures/UC_variancevsDAfreq_depth_facets_IC.jpg"),width=3.5, height=4)
