#Figures for BVR FLARE ms
#09 Sep 2022 HLW

#load libraries
pacman::p_load(dplyr,readr,ggplot2, FSA, AnalystHelper, rcompanion, rstatix, ggpubr, stringr)


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
all_DA_forecasts$phen <- ifelse(all_DA_forecasts$datetime <= as.POSIXct(strat_date) & 
                                  all_DA_forecasts$datetime >="2021-03-13","Stratified", "Mixed")

#remove n=6 days with ice-cover 
all_DA_forecasts <- all_DA_forecasts[!(as.Date(all_DA_forecasts$datetime) %in% c(as.Date("2021-01-10"), as.Date("2021-01-11"),as.Date("2021-01-30"),
                                                                                 as.Date("2021-02-13"),as.Date("2021-02-14"),as.Date("2021-02-15"))),]

#drop 11m completely because some rows were NA when water level was low
all_DA_forecasts <- all_DA_forecasts[!(all_DA_forecasts$depth==11),]

#change model_id to be all uppercase
all_DA_forecasts$model_id <- str_to_title(all_DA_forecasts$model_id)

#only keep 2021 data
all_DA_forecasts <- all_DA_forecasts[all_DA_forecasts$datetime<="2021-12-31",]

#------------------------------------------------------------------------------#
#calculate forecast skill metrics

#forecast skill for each depth and horizon
forecast_skill_depth_horizon <-  plyr::ddply(all_DA_forecasts, c("depth","horizon", "phen", "model_id"), function(x) {
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

#another df for horizon forecast skill aggregated across depths
forecast_horizon_avg <- plyr::ddply(all_DA_forecasts, c("horizon", "model_id", "phen"), function(x) {
  data.frame(
    RMSE = sqrt(mean((x$mean - x$observation)^2, na.rm = TRUE)),
    MAE = mean(abs(x$mean - x$observation), na.rm = TRUE),
    pbias = 100 * (sum(x$mean - x$observation, na.rm = TRUE) / sum(x$observation, na.rm = TRUE)),
    CRPS = verification::crps(x$observation, as.matrix(x[, c(7,9)]))$CRPS,
    variance = (mean(x$sd))^2
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

#order DA frequencies
forecast_horizon_avg$model_id <- factor(forecast_horizon_avg$model_id, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))


#dataframes looking at RMSE across depths, horizons, and periods separately
#forecast_skill_depth <-  plyr::ddply(all_DA_forecasts, c("depth", "model_id"), function(x) {
#  data.frame(
#    RMSE = sqrt(mean((x$mean - x$observation)^2, na.rm = TRUE)),
#    MAE = mean(abs(x$mean - x$observation), na.rm = TRUE),
#    pbias = 100 * (sum(x$mean - x$observation, na.rm = TRUE) / sum(x$observation, na.rm = TRUE)),
#    CRPS = verification::crps(x$observation, as.matrix(x[, c(7,9)]))$CRPS,
#    variance = (mean(x$sd))^2
#  )
#}, .progress = plyr::progress_text(), .parallel = FALSE) 
#
##order DA frequencies
#forecast_skill_depth$model_id <- factor(forecast_skill_depth$model_id, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))

#df with averaged forecast skill for all days (group by horizon, DA, and phen)
#forecast_skill_horizon <- plyr::ddply(all_DA_forecasts, c("horizon", "model_id"), function(x) {
#  data.frame(
#    RMSE = sqrt(mean((x$mean - x$observation)^2, na.rm = TRUE)),
#    MAE = mean(abs(x$mean - x$observation), na.rm = TRUE),
#    pbias = 100 * (sum(x$mean - x$observation, na.rm = TRUE) / sum(x$observation, na.rm = TRUE)),
#    CRPS = verification::crps(x$observation, as.matrix(x[, c(7,9)]))$CRPS,
#    variance = (mean(x$sd))^2
#  )
#}, .progress = plyr::progress_text(), .parallel = FALSE) 
#
##order DA frequencies
#forecast_skill_horizon$model_id <- factor(forecast_skill_horizon$model_id, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))

#df with averaged forecast skill for all days (group by horizon, DA, and phen)
#forecast_skill_phen <- plyr::ddply(all_DA_forecasts, c("phen", "model_id"), function(x) {
#  data.frame(
#    RMSE = sqrt(mean((x$mean - x$observation)^2, na.rm = TRUE)),
#    MAE = mean(abs(x$mean - x$observation), na.rm = TRUE),
#    pbias = 100 * (sum(x$mean - x$observation, na.rm = TRUE) / sum(x$observation, na.rm = TRUE)),
#    CRPS = verification::crps(x$observation, as.matrix(x[, c(7,9)]))$CRPS,
#    variance = (mean(x$sd))^2
#  )
#}, .progress = plyr::progress_text(), .parallel = FALSE) 
#
##order DA frequencies
#forecast_skill_phen$model_id <- factor(forecast_skill_phen$model_id, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))

#------------------------------------------------------------------------------#
#FIGURES
cb_friendly <- c("#117733", "#332288","#AA4499", "#44AA99", "#999933", "#661100")
cb_friendly_2 <- c("#8C510A", "#BF812D", "#C7EAE5", "#35978F") #"#DFC27D", "#DEDEDE",

#FIGURE 6: DA frequency boxplots for 1, 5, and 9m and 1, 7, and 35-day horizons

#creating smaller dataset for kw test w/ 1,5,9m and 1,7,35 days
kw_horizons <- forecast_skill_depth_horizon #[forecast_skill_depth_horizon$depth %in% c(1,5,9) & forecast_skill_depth_horizon$horizon %in% c(1,7,35),]

#kruskal wallis and dunn tests for 1m stratified rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==1])
dunn_strat_1m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1] ~ 
                            kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==1])
rslt_strat_1m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_1m$res, threshold = 0.05)$Letter[c(1,4,2,3)])
rslt_strat_1m <- c("b","ab","a","ab")
median_strat_1m_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id=="Daily"])
median_strat_1m_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id=="Weekly"])
median_strat_1m_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id=="Fortnightly"])
median_strat_1m_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 5m stratified rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==5]))
dunn_strat_5m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5] ~ 
                            as.factor(kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==5]))
rslt_strat_5m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_5m$res, threshold = 0.05)$Letter[c(1,4,2,3)])
rslt_strat_5m <- c("b","a","ab","c")
median_strat_5m_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5  & kw_horizons$model_id=="Daily"])
median_strat_5m_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5  & kw_horizons$model_id=="Weekly"])
median_strat_5m_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5  & kw_horizons$model_id=="Fortnightly"])
median_strat_5m_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==5  & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 9m stratified rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==9]))
dunn_strat_9m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9] ~ 
                            as.factor(kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==9]))
rslt_strat_9m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_9m$res, threshold = 0.05)$Letter[c(1,4,2,3)])
rslt_strat_9m <- c("a","c","b","d")
median_strat_9m_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id=="Daily"])
median_strat_9m_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id=="Weekly"])
median_strat_9m_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id=="Fortnightly"])
median_strat_9m_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id=="Monthly"])
                           
#kruskal wallis and dunn tests for 1m mixed rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==1]))
dunn_mix_1m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1] ~ 
                          as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==1]))
rslt_mix_1m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_1m$res, threshold = 0.05)$Letter[c(1,4,2,3)])
rslt_mix_1m <- c("b","a","ab","ab")
median_mix_1m_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id=="Daily"])
median_mix_1m_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id=="Weekly"])
median_mix_1m_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id=="Fortnightly"])
median_mix_1m_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 5m mixed rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==5]))
dunn_mix_5m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5] ~ 
                          as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==5]))
rslt_mix_5m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_5m$res, threshold = 0.05)$Letter[c(1,4,2,3)])
rslt_mix_5m <- c("b","a","ab","ab")
median_mix_5m_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id=="Daily"])
median_mix_5m_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id=="Weekly"])
median_mix_5m_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id=="Fortnightly"])
median_mix_5m_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 9m mixed rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==9]))
dunn_mix_9m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9] ~ 
                          as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==9]))
rslt_mix_9m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_9m$res, threshold = 0.05)$Letter[c(1,4,2,3)])
rslt_mix_9m <- c("b","a","ab","ab")
median_mix_9m_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id=="Daily"])
median_mix_9m_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id=="Weekly"])
median_mix_9m_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id=="Fortnightly"])
median_mix_9m_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id=="Monthly"])

#Now aggregate across depths
#kruskal wallis and dunn tests for 1d stratified rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1])
dunn_strat_1d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1] ~ 
                            kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1])
rslt_strat_1d=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_1d$res, threshold = 0.05)$Letter[c(1,4,2,3)])
median_strat_1d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==1 & kw_horizons$model_id=="Daily"])
median_strat_1d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==1 & kw_horizons$model_id=="Weekly"])
median_strat_1d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==1 & kw_horizons$model_id=="Fortnightly"])
median_strat_1d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==1 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 7d stratified rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7]))
dunn_strat_7d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7] ~ 
                            as.factor(kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7]))
rslt_strat_7d=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_7d$res, threshold = 0.05)$Letter[c(1,4,2,3)])
median_strat_7d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==7 & kw_horizons$model_id=="Daily"])
median_strat_7d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==7 & kw_horizons$model_id=="Weekly"])
median_strat_7d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==7 & kw_horizons$model_id=="Fortnightly"])
median_strat_7d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==7 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 35d stratified rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35]))
dunn_strat_35d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35] ~ 
                             as.factor(kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35]))
rslt_strat_35d=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_35d$res, threshold = 0.05)$Letter[c(1,4,2,3)])
median_strat_35d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==35 & kw_horizons$model_id=="Daily"])
median_strat_35d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==35 & kw_horizons$model_id=="Weekly"])
median_strat_35d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==35 & kw_horizons$model_id=="Fortnightly"])
median_strat_35d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==35 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 1d mixed rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1]))
dunn_mix_1d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1] ~ 
                          as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1]))
rslt_mix_1d=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_1d$res, threshold = 0.05)$Letter[c(1,4,2,3)])
median_mix_1d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==1 & kw_horizons$model_id=="Daily"])
median_mix_1d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==1 & kw_horizons$model_id=="Weekly"])
median_mix_1d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==1 & kw_horizons$model_id=="Fortnightly"])
median_mix_1d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==1 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 7d mixed rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7]))
dunn_mix_7d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7] ~ 
                          as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7]))
rslt_mix_7d=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_7d$res, threshold = 0.05)$Letter[c(1,4,2,3)])
rslt_mix_7d <- c("ab","a","bc","c")
median_mix_7d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==7 & kw_horizons$model_id=="Daily"])
median_mix_7d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==7 & kw_horizons$model_id=="Weekly"])
median_mix_7d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==7 & kw_horizons$model_id=="Fortnightly"])
median_mix_7d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==7 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 35d mixed rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35]))
dunn_mix_35d <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35] ~ 
                           as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35]))
rslt_mix_35d=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_35d$res, threshold = 0.05)$Letter[c(1,4,2,3)])
rslt_mix_35d <- c("c","a","b","ab")
median_mix_35d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==35 & kw_horizons$model_id=="Daily"])
median_mix_35d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==35 & kw_horizons$model_id=="Weekly"])
median_mix_35d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==35 & kw_horizons$model_id=="Fortnightly"])
median_mix_35d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==35 & kw_horizons$model_id=="Monthly"])

#create table to export mixed and stratified p-vals across depths (aggregated over horizon)
pvals_horizon_aggregated <- data.frame("Depth_m"= c(rep(1,6),rep(5,6),rep(9,6)), 
                                       "Comparison" = c(dunn_mix_1m$res$Comparison,dunn_mix_5m$res$Comparison,dunn_mix_9m$res$Comparison),
                                       "Mixed_pvalue" = c(dunn_mix_1m$res$P.adj,dunn_mix_5m$res$P.adj,dunn_mix_9m$res$P.adj),
                                       "Stratified_pvalue" = c(dunn_strat_1m$res$P.adj,dunn_strat_5m$res$P.adj,dunn_strat_9m$res$P.adj))

#write.csv(pvals_horizon_aggregated,file.path(lake_directory,"analysis/data/pvals_horizon_aggregatred.csv"),row.names = FALSE)

pvals_depth_aggregated <- data.frame("Horizon_days"= c(rep(1,6),rep(7,6),rep(35,6)), 
                                     "Comparison" = c(dunn_mix_1d$res$Comparison,dunn_mix_7d$res$Comparison,dunn_mix_35d$res$Comparison),
                                     "Mixed_pvalue" = c(dunn_mix_1d$res$P.adj,dunn_mix_7d$res$P.adj,dunn_mix_35d$res$P.adj),
                                     "Stratified_pvalue" = c(dunn_strat_1d$res$P.adj,dunn_strat_7d$res$P.adj,dunn_strat_35d$res$P.adj))

#write.csv(pvals_depth_aggregated,file.path(lake_directory,"analysis/data/pvals_depth_aggregatred.csv"),row.names = FALSE)


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
#write.csv(median_RMSE_horizon,file.path(lake_directory,"analysis/data/median_RMSE_depth_horizon_DA.csv"),row.names = FALSE)

#add binned RMSE for new fig 
median_RMSE_horizon$RMSE_bins <- ifelse(median_RMSE_horizon$RMSE_C < 1, 1, # < 1
                                        ifelse(median_RMSE_horizon$RMSE_C >= 1 & median_RMSE_horizon$RMSE_C < 1.5, 1.5, # 1 - 1.5
                                        ifelse(median_RMSE_horizon$RMSE_C >= 1.5 & median_RMSE_horizon$RMSE_C < 2, 2, # 1.5 - 2
                                               ifelse(median_RMSE_horizon$RMSE_C >= 2 & median_RMSE_horizon$RMSE_C < 2.5, 2.5, # 2 - 2.5
                                                      ifelse(median_RMSE_horizon$RMSE_C >= 2.5, 3, 0))))) # 2.5 - 3 

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$TempDynamics=="Mixed"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$TempDynamics=="Stratified"])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$TempDynamics=="Mixed" & median_RMSE_horizon$model_id=="Daily"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$TempDynamics=="Mixed" & median_RMSE_horizon$model_id=="Weekly"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$TempDynamics=="Mixed" & median_RMSE_horizon$model_id=="Fortnightly"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$TempDynamics=="Mixed" & median_RMSE_horizon$model_id=="Monthly"])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Daily"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Weekly"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Fortnightly"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Monthly"])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==1])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==7])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==35])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==1])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==5])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==9])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==1 & median_RMSE_horizon$TempDynamics=="Mixed"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==5 & median_RMSE_horizon$TempDynamics=="Mixed"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==9 & median_RMSE_horizon$TempDynamics=="Mixed"])

range(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==1 & median_RMSE_horizon$TempDynamics=="Stratified"])
range(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==5 & median_RMSE_horizon$TempDynamics=="Stratified"])
range(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==9 & median_RMSE_horizon$TempDynamics=="Stratified"])

median(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==1 & median_RMSE_horizon$TempDynamics=="Mixed"])
median(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==7 & median_RMSE_horizon$TempDynamics=="Mixed"])
median(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==35 & median_RMSE_horizon$TempDynamics=="Mixed"])

median(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==1 & median_RMSE_horizon$TempDynamics=="Stratified"])
median(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==7 & median_RMSE_horizon$TempDynamics=="Stratified"])
median(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==35 & median_RMSE_horizon$TempDynamics=="Stratified"])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==1 & median_RMSE_horizon$TempDynamics=="Mixed" & median_RMSE_horizon$model_id=="Monthly"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==7 & median_RMSE_horizon$TempDynamics=="Mixed" & median_RMSE_horizon$model_id=="Fortnightly"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==35 & median_RMSE_horizon$TempDynamics=="Mixed" & median_RMSE_horizon$model_id=="Weekly"])

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


mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==1 & median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Monthly"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==5 & median_RMSE_horizon$TempDynamics=="Stratified"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Depth_m==9 & median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Monthly"])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==1 & median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Monthly"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==7 & median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Monthly"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==35 & median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Monthly"])

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


# 1-day ahead mixed and strat and 7-day ahead mixed are sig
test <- kw_horizons %>% 
  group_by(horizon,phen) %>% 
  kruskal_test(RMSE ~ model_id)


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#### FIGURE 6: 1,5,9m depth facets for each DA frequency for 1,7,and 35-day ahead forecasts  ####
letters <- data.frame("depth"= c(rep(1,8),rep(5,8),rep(9,8)),
                      "horizon" = rep(c(1,7,35),8),
                      "model_id" = rep(c("Daily","Weekly","Fortnightly","Monthly"),6),
                      "x" = rep(c(1,2,3,4),6),
                      "phen" = rep(c(rep("Mixed",4),rep("Stratified",4)),3),
                      "letter" = tolower(c(rslt_mix_1m, rslt_strat_1m, rslt_mix_5m, rslt_strat_5m,
                                           rslt_mix_9m,rslt_strat_9m)),
                      "max.RMSE" = c(rep(2.6,8),rep(2.5,8),rep(1.6,8)))

#rename depth facets
depths <- c("1m","5m","9m")
names(depths) <- c("1","5","9")

#order factor levels
kw_horizons$model_id <- factor(kw_horizons$model_id, levels = c("Daily", "Weekly", "Fortnightly", "Monthly"))

  ggplot(subset(kw_horizons, horizon %in% c(1,7,35) & depth %in% c(1,5,9)),
         aes(model_id, RMSE, fill=as.factor(horizon))) +  ylab("RMSE") + xlab("")+
  geom_bar(stat="identity",position="dodge") + theme_bw() + guides(fill=guide_legend(title="Horizon")) +
  geom_hline(yintercept=2, linetype='dashed', col = 'black', linewidth=0.3)+ 
  scale_fill_manual(values=c("#81A665","#E0CB48","#D08151"),labels = c("1day", "7day", "35day")) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = c(0.75,0.31),
        legend.background = element_blank(),legend.direction = "horizontal", panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0.05,-0.2,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 6),legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  geom_text(data=letters,aes(x=model_id,y=0.2+max.RMSE,label=letter),hjust=0.1,vjust = -0.1, size=2.5) +
  facet_grid(depth~phen, scales="free_y",labeller = labeller(depth = depths)) 
ggsave(file.path(lake_directory,"analysis/figures/RMSEvsDAfreq_depth_facets_fig5.jpg"),width=3.5, height=4)

#NEW FIGURE 5!!
#now make a plot for RMSE vs all horizons
ggplot(subset(forecast_skill_depth_horizon, depth %in% c(1,5,9)) ,aes(horizon, RMSE, color=as.factor(model_id))) +  ylab("RMSE") + xlab("Horizon (days)")+
   geom_line() + theme_bw() + guides(fill=guide_legend(title="")) + geom_hline(yintercept = 2, linetype="dotted") +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = c(0.89,0.75),
        legend.background = element_blank(),legend.direction = "vertical", panel.grid.minor = element_blank(),
        legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_blank(),legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  facet_grid(depth~phen, scales="free",labeller = labeller(depth = depths)) + scale_color_manual(values=cb_friendly_2) 
ggsave(file.path(lake_directory,"analysis/figures/RMSEvshorizon_depth_facets.jpg"),width=3.5, height=4)

forecast_skill_depth_horizon[forecast_skill_depth_horizon$phen=="Mixed" #& forecast_skill_depth_horizon$horizon ==5
                             &forecast_skill_depth_horizon$model_id=="Monthly" & forecast_skill_depth_horizon$depth==9,]

forecast_skill_depth_horizon[forecast_skill_depth_horizon$phen=="Stratified" #& forecast_skill_depth_horizon$horizon ==5
                             &forecast_skill_depth_horizon$model_id=="Monthly" & forecast_skill_depth_horizon$depth==9,]


#---------------------------------------------------------------------------------#
#### SI figure: mixed and stratified RMSE tileplots aggregated across all depths  ####

#round rmse to nearest 0.5 for tile plot below
forecast_horizon_avg$RMSE_bins <- plyr::round_any(forecast_horizon_avg$RMSE,0.5) 

#figure for horizon vs frequency to compare forecast skill
ggplot(forecast_horizon_avg, aes(model_id, horizon, fill=RMSE_bins)) + 
  ylab("Horizon (days)") + theme_bw() + facet_wrap(~phen) + geom_tile(width=0.8) + xlab("") +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.spacing=unit(0, "cm"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  guides(fill=guide_legend(title="RMSE")) +  scale_fill_gradientn(labels=c(
    "< 0.5","0.5 - 1.0","1.0 - 1.5","1.5 - 2.0","2.0 - 2.5","2.5 - 3.0"),colors = hcl.colors(4, "BuPu")) 
ggsave(file.path(lake_directory,"analysis/figures/HorizonvsDA_tileplot_fig6.jpg"),width=4, height=3.5)

#median mixed vs stratified period rmse across different horizons
median(forecast_horizon_avg$RMSE_bins[forecast_horizon_avg$phen=="Mixed" & forecast_horizon_avg$horizon==5 & forecast_horizon_avg$model_id=="Weekly"])
median(forecast_horizon_avg$RMSE_bins[forecast_horizon_avg$phen=="Stratified" & forecast_horizon_avg$horizon==5 & forecast_horizon_avg$model_id=="Weekly"])

#------------------------------------------------------------------------------------------------#
# Data assimilation figure 4
DA <- all_DA_forecasts

#pull out 2 horizons for fig 4 (mixed vs stratified)
DA_sub <- DA[(DA$datetime >="2021-01-01" & DA$datetime <="2021-01-31") | (DA$datetime >="2021-06-24" & DA$datetime <="2021-07-25"),]

#summary df to average the forecasts for each DA freq, horizon, depth, and date
DA_sub_final <- plyr::ddply(DA_sub, c("depth","horizon","datetime", "model_id"), function(x) {
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
ggsave(file.path(lake_directory,"analysis/figures/AssimilationVSdatafreq_fig4.jpg"), width=3.5, height=4) 

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

#-------------------------------------------------------------------------------#
# Fig 3 - phenology plot w/ temp at different depths 
ggplot(sub) +   theme_bw() +
  geom_rect(data = phen, aes(fill = "Mixed"), xmin=-Inf ,xmax = as.Date("2021-03-12"), ymin = -Inf, ymax = Inf, inherit.aes = FALSE) + 
  geom_rect(data = phen, aes(fill = "Stratified"), xmin=as.Date("2021-03-13") ,xmax = as.Date("2021-11-07"), ymin = -Inf, ymax = Inf, inherit.aes = FALSE)+
  geom_rect(data = phen, aes(fill = "Mixed"), xmin=as.Date("2021-11-08") ,xmax = Inf, ymin = -Inf, ymax = Inf, inherit.aes = FALSE) +
  geom_line(aes(Date, as.numeric(observation), color = factor(depth)), size=0.4) + ylab(expression("Temperature ("*~degree*C*")")) + xlab("") +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = c(0.11,0.59), legend.background = element_blank(),
        legend.key = element_blank(), legend.key.height = unit(0.3,"cm"), legend.key.width = unit(0.4,"cm"), legend.spacing.y = unit(0.01,"cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_fill_manual('', values = c('gray','white')) +
  guides(color = guide_legend("Depth (m)"), fill= guide_legend(order = 1, override.aes= list(color="black")))
ggsave(file.path(lake_directory,"analysis/figures/2021_watertemp_mixedVstratified.jpg"))

#-------------------------------------------------------------------------------#
# Fig 7 - fig to compare forecast skill across horizons and depths 

#test normality √ (but variances are not homogenous so stick w/ KW)
ggpubr::ggqqplot(kw_horizons$RMSE)

#create new df with all letters 
DA_horizon_letters <- data.frame("Horizon" = rep(c(rep(1,4),rep(7,4), rep(35,4)),2),
                                 "TempDynamics" = c(rep("Mixed",12),rep("Stratified",12)),
                                 "model_id" = c(rep(c("Daily","Weekly","Fortnightly","Monthly"),6)),
                                 "RMAE" = tolower(c(rslt_mix_1d,rslt_mix_7d,rslt_mix_35d,
                                                      rslt_strat_1d,rslt_strat_7d,rslt_strat_35d)))
                                                    
DA_depth_letters <- data.frame("Depth" = rep(c(rep(1,4),rep(5,4), rep(9,4)),2),
                                 "TempDynamics" = c(rep("Mixed",12),rep("Stratified",12)),
                                 "model_id" = c(rep(c("Daily","Weekly","Fortnightly","Monthly"),6)),
                                 "letter" = tolower(c(rslt_mix_1m,rslt_mix_5m,rslt_mix_9m,
                                                      rslt_strat_1m,rslt_strat_5m,rslt_strat_9m)))



#change factor order
DA_horizon_letters$model_id <- factor(DA_horizon_letters$model_id, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))
DA_depth_letters$model_id <- factor(DA_depth_letters$model_id, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))

#change horizon order
median_RMSE_horizon$Horizon_days <- factor(median_RMSE_horizon$Horizon_days, levels=c(1,7,35))
median_RMSE_horizon$Depth_m <- factor(median_RMSE_horizon$Depth_m, levels=c(9,5,1))

#now make a figure... 
horiz <- ggplot(median_RMSE_horizon, aes(model_id, as.factor(Horizon_days), fill=as.factor(RMSE_bins))) + 
  ylab("Horizon (days)") + xlab("") + theme_bw() + geom_tile(width=0.8) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = "right",
        legend.background = element_blank(),legend.direction = "vertical", panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0.05,0.2,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 6),legend.text  = element_text(size = 8), panel.spacing=unit(0, "cm"), legend.margin=margin(0,0,0,0),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(size=6)) +
  facet_grid(~TempDynamics, scales="free_y",labeller = labeller(Depth_m = depths)) +
  guides(fill=guide_legend(title="RMSE")) + scale_fill_brewer(type="seq", 
    palette=7, direction = -1, labels = c("< 1", "1 - 1.5", "1.5 - 2", "2 - 2.5", "> 2.5"))

depth <- ggplot(median_RMSE_horizon, aes(model_id, as.factor(Depth_m), fill=as.factor(RMSE_bins))) + 
  ylab("Depth (m)") + xlab("") + theme_bw() + geom_tile(width=0.8) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = "right",
        legend.background = element_blank(),legend.direction = "vertical", panel.grid.minor = element_blank(),
        plot.margin = unit(c(-0.6,0.05,-0.2,0.1), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 6),legend.text  = element_text(size = 8), panel.spacing=unit(0, "cm"), legend.margin=margin(0,0,0,0),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  facet_grid(~TempDynamics, scales="free_y",labeller = labeller(Depth_m = depths)) +
  guides(fill=guide_legend(title="RMSE")) + scale_fill_brewer(type="seq", 
    palette=7, direction = -1, labels = c("< 1", "1 - 1.5", "1.5 - 2", "2 - 2.5", "> 2.5"))

#combine depth and horizon plots
ggarrange(horiz, depth, 
          ncol = 1, nrow = 2,
          common.legend = TRUE,
          legend="right")

ggsave(file.path(lake_directory,"analysis/figures/RMSE_bins_horizon_depth_facets_fig7.jpg"),width=3.5, height=4)
