#Figures for BVR FLARE ms
#09 Sep 2022 HLW

#load libraries
pacman::p_load(dplyr,readr,ggplot2, FSA, AnalystHelper, rcompanion, rstatix, ggpubr)

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
score_dir <- arrow::SubTreeFileSystem$create(file.path(lake_directory,"scores/da_study"))
all_DA_forecasts <- arrow::open_dataset(score_dir) |> collect() |>   
  filter(!is.na(observation), variable == "temperature",horizon > 0)

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
all_DA_forecasts$phen <- ifelse(all_DA_forecasts$datetime <= as.Date(strat_date) & 
                                  all_DA_forecasts$datetime >="2021-03-13","Stratified", "Mixed")

#remove n=6 days with ice-cover 
all_DA_forecasts <- all_DA_forecasts[!(as.Date(all_DA_forecasts$datetime) %in% c(as.Date("2021-01-10"), as.Date("2021-01-11"),as.Date("2021-01-30"),
                                    as.Date("2021-02-13"),as.Date("2021-02-14"),as.Date("2021-02-15"))),]

#change model_id to be all uppercase
all_DA_forecasts$model_id <- str_to_title(all_DA_forecasts$model_id)
#------------------------------------------------------------------------------#
#calculate forecast skill metrics

#forecast skill for each depth and horizon
forecast_skill_depth_horizon <-  plyr::ddply(all_DA_forecasts, c("depth", "datetime","horizon", "model_id"), function(x) {
  data.frame(
    RMSE = sqrt(mean((x$mean - x$observation)^2, na.rm = TRUE)),
    MAE = mean(abs(x$mean - x$observation), na.rm = TRUE),
    pbias = 100 * (sum(x$mean - x$observation, na.rm = TRUE) / sum(x$observation, na.rm = TRUE)),
    CRPS = verification::crps(x$observation, as.matrix(x[, c(7,9)]))$CRPS
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

##add in mixed/stratified period
forecast_skill_depth_horizon$phen <- ifelse(as.Date(forecast_skill_depth_horizon$datetime) <= as.Date(strat_date) & 
                                              as.Date(forecast_skill_depth_horizon$datetime) >="2021-03-13","Stratified", "Mixed")

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

#------------------------------------------------------------------------------#
#FIGURES
cb_friendly <- c("#117733", "#332288","#AA4499", "#44AA99", "#999933", "#661100")
cb_friendly_2 <- c("#8C510A", "#BF812D", "#C7EAE5", "#35978F") #"#DFC27D", "#DEDEDE",

#FIGURE 4: DA frequency boxplots for 1, 5, and 9m and 1, 7, and 35-day horizons

#creating smaller dataset for kw test w/ 1,5,9m and 1,7,35 days
kw_horizons <- forecast_skill_depth_horizon[forecast_skill_depth_horizon$depth %in% c(1,5,9) & forecast_skill_depth_horizon$horizon %in% c(1,7,35),]

#kruskal wallis and dunn tests for 1m stratified rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==1])
dunn_strat_1m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1] ~ 
                            kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==1])
rslt_strat_1m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_1m$res, threshold = 0.05)$Letter[c(1,4,2,3)])
rslt_strat_1m <- c("b","a","a","a") #manually creating letters because I seriously can't figure out how to get a = lowest RMSE
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

median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id=="Daily"])




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

median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id=="Daily"])




#kruskal wallis and dunn tests for 9m stratified rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==9]))
dunn_strat_9m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9] ~ 
                            as.factor(kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==9]))
rslt_strat_9m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_9m$res, threshold = 0.05)$Letter[c(1,4,2,3)])
rslt_strat_9m <- c("c","a","ab","b")
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

median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id=="Daily"])



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

median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id=="Daily"])



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

median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id=="Daily"])



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

median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id=="Daily"])



#create table to export mixed and stratified p-vals across depths (aggregated over horizon)
pvals_horizon_aggregated <- data.frame("Depth_m"= c(rep(1,6),rep(5,6),rep(9,6)), 
                                                  "Comparison" = c(dunn_mix_1m$res$Comparison,dunn_mix_5m$res$Comparison,dunn_mix_9m$res$Comparison),
                                                  #"Z" = c(dunn_mix_1m$res$Z,dunn_mix_5m$res$Z,dunn_mix_9m$res$Z),
                                                  "Mixed_pvalue" = c(dunn_mix_1m$res$P.adj,dunn_mix_5m$res$P.adj,dunn_mix_9m$res$P.adj),
                                                  "Stratified_pvalue" = c(dunn_strat_1m$res$P.adj,dunn_strat_5m$res$P.adj,dunn_strat_9m$res$P.adj))

#write.csv(pvals_horizon_aggregated,file.path(lake_directory,"analysis/data/pvals_horizon_aggregatred.csv"),row.names = FALSE)

#median RMSE table
median_RMSE_horizon <- data.frame("Depth_m" = c(rep(1,24),rep(5,24),rep(9,24)),
                                  "DA" = rep(c(rep("Daily",3), rep("Weekly",3),rep("Fortnightly",3),rep("Monthly",3)),6),
                                  "Horizon_days" = rep(c(1,7,35),24),
                                  "TempDynamics" = rep(c(rep("Mixed",12),rep("Stratified",12)),3),
                                  "RMSE_C" = c(median_mix_1m_daily,median_mix_1m_weekly,median_mix_1m_fortnightly,median_mix_1m_monthly,
                                             median_strat_1m_daily,median_strat_1m_weekly,median_strat_1m_fortnightly,median_strat_1m_monthly,
                                             median_mix_5m_daily,median_mix_5m_weekly,median_mix_5m_fortnightly,median_mix_5m_monthly,
                                             median_strat_5m_daily,median_strat_5m_weekly,median_strat_5m_fortnightly,median_strat_5m_monthly,
                                             median_mix_9m_daily,median_mix_9m_weekly,median_mix_9m_fortnightly,median_mix_9m_monthly,
                                             median_strat_9m_daily,median_strat_9m_weekly,median_strat_9m_fortnightly,median_strat_9m_monthly))
#write.csv(median_RMSE_horizon,file.path(lake_directory,"analysis/data/median_RMSE_depth_horizon_DA.csv"),row.names = FALSE)

#max(median_RMSE_horizon$RMSE[median_RMSE_horizon$Depth_m==1 & median_RMSE_horizon$TempDynamics=="Stratified"]) - 
#  min(median_RMSE_horizon$RMSE[median_RMSE_horizon$depth==1 & median_RMSE_horizon$TempDynamics=="Stratified"])
#max(median_RMSE_horizon$RMSE[median_RMSE_horizon$depth==5 & median_RMSE_horizon$TempDynamics=="Stratified"]) -
#  min(median_RMSE_horizon$RMSE[median_RMSE_horizon$depth==5 & median_RMSE_horizon$TempDynamics=="Stratified"])
#max(median_RMSE_horizon$RMSE[median_RMSE_horizon$depth==9 & median_RMSE_horizon$TempDynamics=="Stratified"]) -
#  min(median_RMSE_horizon$RMSE[median_RMSE_horizon$depth==9 & median_RMSE_horizon$TempDynamics=="Stratified"])

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
rslt_strat_weekly_1m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_weekly_1m$res, threshold = 0.05)$Letter)
rslt_strat_weekly_1m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==1], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==7], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly"]))
dunn_strat_fortnightly_1m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly"] ~ 
                                   as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly"]))
rslt_strat_fortnightly_1m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_fortnightly_1m$res, threshold = 0.05)$Letter)
rslt_strat_fortnightly_1m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==1], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==7], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly"]))
dunn_strat_monthly_1m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly"] ~ 
                                        as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly"]))
rslt_strat_monthly_1m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_monthly_1m$res, threshold = 0.05)$Letter)
rslt_strat_monthly_1m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly" & kw_horizons$horizon ==1], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly" & kw_horizons$horizon ==7], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly" & kw_horizons$horizon ==35], probs=.75))


#kruskal wallis and dunn tests for 1-day ahead stratified rmse across depth and DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Daily"] ~ 
               as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Daily"]))
dunn_strat_daily_1day <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Daily"] ~ 
                                  as.factor(kw_horizons$depth[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1 & kw_horizons$model_id == "Daily"]))
rslt_strat_daily_1day=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_daily_1day$res, threshold = 0.05)$Letter[c(1,3,2)])
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
rslt_strat_daily_5m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_daily_5m$res, threshold = 0.05)$Letter)
rslt_strat_daily_5m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==1], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==7], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly"]))
dunn_strat_weekly_5m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly"] ~ 
                                   as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly"]))
rslt_strat_weekly_5m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_weekly_5m$res, threshold = 0.05)$Letter)
rslt_strat_weekly_5m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==1], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==7], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly"]))
dunn_strat_fortnightly_5m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly"] ~ 
                                        as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly"]))
rslt_strat_fortnightly_5m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_fortnightly_5m$res, threshold = 0.05)$Letter)
rslt_strat_fortnightly_5m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==1], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==7], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Monthly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Monthly"]))
dunn_strat_monthly_5m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Monthly"] ~ 
                                    as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==5 & kw_horizons$model_id == "Monthly"]))
rslt_strat_monthly_5m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_monthly_5m$res, threshold = 0.05)$Letter)
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
rslt_strat_daily_9m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_daily_9m$res, threshold = 0.05)$Letter)
rslt_strat_daily_9m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==1], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==7], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly"]))
dunn_strat_weekly_9m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly"] ~ 
                                   as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly"]))
rslt_strat_weekly_9m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_weekly_9m$res, threshold = 0.05)$Letter)
rslt_strat_weekly_9m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==1], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==7], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly"]))
dunn_strat_fortnightly_9m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly"] ~ 
                                        as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly"]))
rslt_strat_fortnightly_9m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_fortnightly_9m$res, threshold = 0.05)$Letter)
rslt_strat_fortnightly_9m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==1], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==7], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Monthly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Monthly"]))
dunn_strat_monthly_9m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Monthly"] ~ 
                                    as.factor(kw_horizons$horizon[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id == "Monthly"]))
rslt_strat_monthly_9m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_monthly_9m$res, threshold = 0.05)$Letter)
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
rslt_mix_daily_1m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_daily_1m$res, threshold = 0.05)$Letter)
rslt_mix_daily_1m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==1], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==7], probs=.75),
                             quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly"]))
dunn_mix_weekly_1m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly"] ~ 
                                   as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly"]))
rslt_mix_weekly_1m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_weekly_1m$res, threshold = 0.05)$Letter)
rslt_mix_weekly_1m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==1], probs=.75),
                           quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==7], probs=.75),
                           quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly"]))
dunn_mix_fortnightly_1m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly"] ~ 
                                        as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly"]))
rslt_mix_fortnightly_1m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_fortnightly_1m$res, threshold = 0.05)$Letter)
rslt_mix_fortnightly_1m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==1], probs=.75),
                           quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==7], probs=.75),
                           quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly"]))
dunn_mix_monthly_1m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly"] ~ 
                                    as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id == "Monthly"]))
rslt_mix_monthly_1m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_monthly_1m$res, threshold = 0.05)$Letter)
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
rslt_mix_daily_5m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_daily_5m$res, threshold = 0.05)$Letter)
rslt_mix_daily_5m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==1], probs=.75),
                           quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==7], probs=.75),
                           quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly"]))
dunn_mix_weekly_5m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly"] ~ 
                                   as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly"]))
rslt_mix_weekly_5m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_weekly_5m$res, threshold = 0.05)$Letter)
rslt_mix_weekly_5m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==1], probs=.75),
                           quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==7], probs=.75),
                           quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly"]))
dunn_mix_fortnightly_5m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly"] ~ 
                                        as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly"]))
rslt_mix_fortnightly_5m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_fortnightly_5m$res, threshold = 0.05)$Letter)
rslt_mix_fortnightly_5m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==1], probs=.75),
                           quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==7], probs=.75),
                           quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$v == "Fortnightly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Monthly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Monthly"]))
dunn_mix_monthly_5m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Monthly"] ~ 
                                    as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id == "Monthly"]))
rslt_mix_monthly_5m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_monthly_5m$res, threshold = 0.05)$Letter)
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
rslt_mix_daily_9m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_daily_9m$res, threshold = 0.05)$Letter)
rslt_mix_daily_9m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==1], probs=.75),
                           quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==7], probs=.75),
                           quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Daily" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly"]))
dunn_mix_weekly_9m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly"] ~ 
                                   as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly"]))
rslt_mix_weekly_9m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_weekly_9m$res, threshold = 0.05)$Letter)
rslt_mix_weekly_9m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==1], probs=.75),
                           quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==7], probs=.75),
                           quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Weekly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly"]))
dunn_mix_fortnightly_9m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly"] ~ 
                                        as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly"]))
rslt_mix_fortnightly_9m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_fortnightly_9m$res, threshold = 0.05)$Letter)
rslt_mix_fortnightly_9m_max <- c(quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==1], probs=.75),
                           quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==7], probs=.75),
                           quantile(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Fortnightly" & kw_horizons$horizon ==35], probs=.75))

kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Monthly"] ~ 
               as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Monthly"]))
dunn_mix_monthly_9m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Monthly"] ~ 
                                    as.factor(kw_horizons$horizon[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id == "Monthly"]))
rslt_mix_monthly_9m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_monthly_9m$res, threshold = 0.05)$Letter)
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
#write.csv(pvals,file.path(lake_directory,"analysis/data/pvals_horizon.csv"),row.names = FALSE)

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
#write.csv(pvals_depths,file.path(lake_directory,"analysis/data/pvals_depth.csv"),row.names = FALSE)


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#### FIGURE 4: 1,5,9m depth facets for each DA frequency for 1,7,and 35-day ahead forecasts  ####

#new df with letters for each horizon
#letters <- data.frame("depth"= c(rep(1,24),rep(5,24),rep(9,24)),
#                      "horizon" = rep(c(1,7,35),24),
#                      "DA" = rep(c(rep("Daily",3),rep("Weekly",3),rep("Fortnightly",3),rep("Monthly",3)),6),
#                      "x" = rep(c(0.8,1.05,1.3,1.8,2.05,2.3,2.8,3.05,3.3,3.8,4.05,4.3),6),
#                      "phen" = rep(c(rep("Mixed",12),rep("Stratified",12)),3),
#                      "letter" = tolower(c(rslt_mix_daily_1m, rslt_mix_weekly_1m, rslt_mix_fortnightly_1m, rslt_mix_monthly_1m,
#                                            rslt_strat_daily_1m, rslt_strat_weekly_1m, rslt_strat_fortnightly_1m, rslt_strat_monthly_1m,
#                                            rslt_mix_daily_5m, rslt_mix_weekly_5m, rslt_mix_fortnightly_5m, rslt_mix_monthly_5m,
#                                            rslt_strat_daily_5m, rslt_strat_weekly_5m, rslt_strat_fortnightly_5m, rslt_strat_monthly_5m,
#                                            rslt_mix_daily_9m, rslt_mix_weekly_9m, rslt_mix_fortnightly_9m, rslt_mix_monthly_9m,
#                                            rslt_strat_daily_9m, rslt_strat_weekly_9m, rslt_strat_fortnightly_9m, rslt_strat_monthly_9m)),
#                      "max.RMSE" = c(rslt_mix_daily_1m_max, rslt_mix_weekly_1m_max, rslt_mix_fortnightly_1m_max, rslt_mix_monthly_1m_max,
#                                       rslt_strat_daily_1m_max, rslt_strat_weekly_1m_max, rslt_strat_fortnightly_1m_max, rslt_strat_monthly_1m_max,
#                                       rslt_mix_daily_5m_max, rslt_mix_weekly_5m_max, rslt_mix_fortnightly_5m_max, rslt_mix_monthly_5m_max,
#                                       rslt_strat_daily_5m_max, rslt_strat_weekly_5m_max, rslt_strat_fortnightly_5m_max, rslt_strat_monthly_5m_max,
#                                       rslt_mix_daily_9m_max, rslt_mix_weekly_9m_max, rslt_mix_fortnightly_9m_max, rslt_mix_monthly_9m_max,
#                                       rslt_strat_daily_9m_max, rslt_strat_weekly_9m_max, rslt_strat_fortnightly_9m_max, rslt_strat_monthly_9m_max))

#change max rmse to 3 for those >3 because removing outliers for plots below
#letters$max.RMSE[which(letters$max.RMSE>=3)] <- 2.6

letters <- data.frame("depth"= c(rep(1,8),rep(5,8),rep(9,8)),
                      "horizon" = rep(c(1,7,35),8),
                      "model_id" = rep(c("Daily","Weekly","Fortnightly","Monthly"),6),
                      "x" = rep(c(1,2,3,4),6),
                      "phen" = rep(c(rep("Mixed",4),rep("Stratified",4)),3),
                      "letter" = tolower(c(rslt_mix_1m, rslt_strat_1m, rslt_mix_5m, rslt_strat_5m,
                                            rslt_mix_9m,rslt_strat_9m)),
                      "max.RMSE" = c(rep(3.5,16),rep(1.6,8)))

#rename depth facets
depths <- c("1m","5m","9m")
names(depths) <- c("1","5","9")

#order factor levels
kw_horizons$model_id <- factor(kw_horizons$model_id, levels = c("Daily", "Weekly", "Fortnightly", "Monthly"))

kw_horizons %>%
  group_by(model_id,depth,horizon,phen) %>%  # do the same calcs for each box
  mutate(value2 = filter_lims(RMSE)) %>%
  ggplot(aes(model_id, value2, fill=as.factor(horizon))) +  ylab("RMSE") + xlab("")+
  geom_boxplot(outlier.shape = NA) + theme_bw() + guides(fill=guide_legend(title="Horizon")) +
  geom_hline(yintercept=2, linetype='dashed', col = 'black', size=0.3)+ 
  scale_fill_manual(values=c("#81A665","#E0CB48","#D08151"),labels = c("1day", "7day", "35day")) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = c(0.75,0.31),
        legend.background = element_blank(),legend.direction = "horizontal", panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0.05,-0.2,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 6),legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  geom_text(data=letters,aes(x=model_id,y=0.2+max.RMSE,label=letters$letter),hjust=0.1,vjust = -0.1, size=2.5) +
  facet_grid(depth~phen, scales="free_y",labeller = labeller(depth = depths)) 
ggsave(file.path(lake_directory,"analysis/figures/RMSEvsDAfreq_depth_facets_fig4.jpg"),width=3.5, height=4)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#### FIGURE 5: mixed and stratified RMSE tileplots aggregated across all depths  ####

#round rmse to nearest 0.5 for tile plot below
forecast_horizon_avg$RMSE_bins <- plyr::round_any(forecast_horizon_avg$RMSE,0.5) 

#figure for horizon vs frequency to compare forecast skill
ggplot(forecast_horizon_avg, aes(model_id, horizon, fill=RMSE_bins)) + 
  ylab("Horizon (days)") + theme_bw() + facet_wrap(~phen) + geom_tile(width=0.8) + xlab("") +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.spacing=unit(0, "cm"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  guides(fill=guide_legend(title="RMSE")) +  scale_fill_gradientn(colors = hcl.colors(4, "BuPu")) 
ggsave(file.path(lake_directory,"analysis/figures/HorizonvsDA_tileplot_fig5.jpg"),width=4, height=3.5)

#median mixed vs stratified period rmse across different horizons
median(forecast_horizon_avg$RMSE_bins[forecast_horizon_avg$phen=="Mixed" & forecast_horizon_avg$horizon<=3])
median(forecast_horizon_avg$RMSE_bins[forecast_horizon_avg$phen=="Stratified" & forecast_horizon_avg$horizon==13 & forecast_horizon_avg$DA=="Daily"])

#------------------------------------------------------------------------------------------------#
# Data assimilation figure 6
DA <- all_DA_forecasts

#pull out 2 horizons for fig 6 (mixed vs stratified)
DA_sub <- DA[(DA$datetime >="2021-01-01" & DA$datetime <="2021-01-31") | (DA$datetime >="2021-06-25" & DA$datetime <="2021-07-25"),]

#summary df to average the forecasts for each DA freq, horizon, depth, forecast_start_day, and date
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
                  color=model_id, fill=model_id, linetype=NA),alpha=0.4,show.legend = F) + geom_line(size=0.5) + theme_bw() + 
  facet_wrap(depth~phen, scales = "free", labeller = labeller(depth=depths,.multi_line = FALSE),ncol = 2) + scale_x_date() +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = c(0.76,0.02),
        legend.background = element_blank(),legend.direction = "horizontal", panel.grid.minor = element_blank(), legend.key=element_rect(fill=NA),
        plot.margin = unit(c(0,0.05,-0.2,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 4.5),legend.text  = element_text(size = 4.5), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  ylab(expression("Temperature ("*~degree*C*")")) + xlab("")  + scale_color_manual(values=rev(cb_friendly_2)) + scale_fill_manual(values=rev(cb_friendly_2)) +
  guides(color = guide_legend(title="",override.aes=list(fill=NA), reverse = TRUE), fill= "none")
ggsave(file.path(lake_directory,"analysis/figures/AssimilationVSdatafreq_fig6.jpg"), width=3.5, height=4) 

#uncertainty range across depths and mixed vs stratified periods

((mean(DA_sub_final$forecast_lower_95[DA_sub_final$phen=="Mixed" & DA_sub_final$model_id=="Monthly"]) -
  mean(DA_sub_final$forecast_upper_95[DA_sub_final$phen=="Mixed" & DA_sub_final$model_id=="Monthly"])) -
  
(mean(DA_sub_final$forecast_lower_95[DA_sub_final$phen=="Mixed" & DA_sub_final$model_id=="Daily"]) -
mean(DA_sub_final$forecast_upper_95[DA_sub_final$phen=="Mixed" & DA_sub_final$model_id=="Daily"]))) /
  
  mean(mean(DA_sub_final$forecast_lower_95[DA_sub_final$phen=="Mixed" & DA_sub_final$model_id=="Monthly"]),
       mean(DA_sub_final$forecast_upper_95[DA_sub_final$phen=="Mixed" & DA_sub_final$model_id=="Monthly"]),
       mean(DA_sub_final$forecast_lower_95[DA_sub_final$phen=="Mixed" & DA_sub_final$model_id=="Daily"]),
         mean(DA_sub_final$forecast_upper_95[DA_sub_final$phen=="Mixed" & DA_sub_final$model_id=="Daily"])) *100


((mean(DA_sub_final$forecast_lower_95[DA_sub_final$phen=="Stratified" & DA_sub_final$model_id=="Monthly"]) -
  mean(DA_sub_final$forecast_upper_95[DA_sub_final$phen=="Stratified" & DA_sub_final$model_id=="Monthly"])) -
     
(mean(DA_sub_final$forecast_lower_95[DA_sub_final$phen=="Stratified" & DA_sub_final$model_id=="Daily"]) -
  mean(DA_sub_final$forecast_upper_95[DA_sub_final$phen=="Stratified" & DA_sub_final$model_id=="Daily"]))) / 
  
  mean(mean(DA_sub_final$forecast_lower_95[DA_sub_final$phen=="Stratified" & DA_sub_final$model_id=="Daily"]),
         mean(DA_sub_final$forecast_upper_95[DA_sub_final$phen=="Stratified" & DA_sub_final$model_id=="Daily"]),
       mean(DA_sub_final$forecast_lower_95[DA_sub_final$phen=="Stratified" & DA_sub_final$model_id=="Monthly"]),
         mean(DA_sub_final$forecast_upper_95[DA_sub_final$phen=="Stratified" & DA_sub_final$model_id=="Monthly"])) *100

#--------------------------------------------------------------------------------#
#Fig 3 - water temp phenology fig
# Read in FLARE observations ----
target_file <- file.path(config$file_path$qaqc_data_directory,
                         "bvre-targets-insitu.csv")
obs <- read.csv(target_file)
wtemp <- obs[obs$variable == "temperature", ] # Subset to temperature

sub <- wtemp[wtemp$date > "2021-01-01" & wtemp$date < "2022-01-02", ] # Subset to target dates - need to change for your experiment
sub$Date <- as.Date(sub$date) # Line
#sub$dens <- rLakeAnalyzer::water.density(sub$value) # density for stratification

# Classifies phenology
phen <- plyr::ddply(sub, "Date", function(x) {
  surf <- x$dens[x$depth == min(x$depth)]
  bott <- x$dens[x$depth == max(x$depth)]
  if(nrow(x) == 1) return(data.frame(status = NA))
  status <- ifelse(abs(surf - bott) < 0.1, "Mixed", "Stratified")
  data.frame(status = status)
})

#round depths up
sub$depth <- floor(sub$depth)


# Fig 3 - phenology plot w/ temp at different depths 
ggplot(sub) +   theme_bw() +
  geom_rect(data = phen, aes(fill = "Mixed"), xmin=-Inf ,xmax = as.Date("2021-03-12"), ymin = -Inf, ymax = Inf, inherit.aes = FALSE) + 
  geom_rect(data = phen, aes(fill = "Stratified"), xmin=as.Date("2021-03-13") ,xmax = as.Date("2021-11-07"), ymin = -Inf, ymax = Inf, inherit.aes = FALSE)+
  geom_rect(data = phen, aes(fill = "Mixed"), xmin=as.Date("2021-11-08") ,xmax = Inf, ymin = -Inf, ymax = Inf, inherit.aes = FALSE) +
  geom_line(aes(Date, as.numeric(observation), color = factor(depth)), size=0.2) + ylab(expression("Temperature ("*~degree*C*")")) + xlab("") +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = c(0.11,0.59), legend.background = element_blank(),
        legend.key = element_blank(), legend.key.height = unit(0.3,"cm"), legend.key.width = unit(0.4,"cm"), legend.spacing.y = unit(0.01,"cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_fill_manual('', values = c('gray','white')) +
  guides(color = guide_legend("Depth (m)"), fill= guide_legend(order = 1, override.aes= list(color="black")))
#ggsave(file.path(lake_directory,"analysis/figures/2021_watertemp_mixedVstratified.jpg"))

