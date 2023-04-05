#Figures for BVR FLARE ms
#09 Sep 2022 HLW

#load libraries
pacman::p_load(dplyr,readr,ggplot2, FSA, AnalystHelper, rcompanion, rstatix, ggpubr, stringr, egg, viridis, padr, ggnewscale, purrr)

#change tag_facet code
tag_facet2 <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                       hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
}

#set wd
lake_directory <- here::here()
setwd(lake_directory)

#read in all forecasts 
score_dir <- arrow::SubTreeFileSystem$create(file.path(lake_directory,"scores/da_study"))
all_DA_forecasts <- arrow::open_dataset(score_dir) |> collect() |>   
  filter(!is.na(observation), variable == "temperature",horizon >= 0, 
         as.Date(reference_datetime) > "2020-12-31") 

#round horizon because Jan 01 2021 has weird decimals
all_DA_forecasts$horizon <- round(all_DA_forecasts$horizon)

#round depths up to nearest m 
#all_DA_forecasts$depth <- ceiling(all_DA_forecasts$depth)

#add a group number so that I can average horizons later on
all_DA_forecasts <- all_DA_forecasts %>% 
  mutate(group = case_when(all_DA_forecasts$horizon <= 5 ~ "1-5",
                           all_DA_forecasts$horizon <=10 & all_DA_forecasts$horizon > 5 ~ "6-10",
                           all_DA_forecasts$horizon <=15 & all_DA_forecasts$horizon > 10 ~ "11-15",
                           all_DA_forecasts$horizon <=20 & all_DA_forecasts$horizon > 15 ~ "16-20",
                           all_DA_forecasts$horizon <=25 & all_DA_forecasts$horizon > 20 ~ "21-25",
                           all_DA_forecasts$horizon <=30 & all_DA_forecasts$horizon > 25 ~ "26-30",
                           all_DA_forecasts$horizon <=35 & all_DA_forecasts$horizon > 30 ~ "31-35"))

strat_date<- "2021-11-07"

#add stratified vs mixed col
all_DA_forecasts$phen <- ifelse(all_DA_forecasts$datetime <= as.POSIXct(strat_date) & 
                                  all_DA_forecasts$datetime >="2021-03-12","Stratified", "Mixed")

#remove n=6 days with ice-cover 
all_DA_forecasts <- all_DA_forecasts[!(as.Date(all_DA_forecasts$datetime) %in% c(as.Date("2021-01-10"), as.Date("2021-01-11"),as.Date("2021-01-30"),
                                                                                 as.Date("2021-02-13"),as.Date("2021-02-14"),as.Date("2021-02-15"))),]

#drop 11m completely because some rows were NA when water level was low
all_DA_forecasts <- all_DA_forecasts[!(all_DA_forecasts$depth==11),]

#change model_id to be all uppercase
all_DA_forecasts$model_id <- str_to_title(all_DA_forecasts$model_id)

#only keep 2021 data
all_DA_forecasts <- all_DA_forecasts[all_DA_forecasts$datetime<="2021-12-31",]

#capitalize first letter in model id
str_sub(all_DA_forecasts$model_id, 1, 1) <- str_sub(all_DA_forecasts$model_id, 1, 1) %>% str_to_upper()

#------------------------------------------------------------------------------#
#quick look to see what water temp ranges during mixed and stratified periods are
range(all_DA_forecasts$observation[all_DA_forecasts$phen=="Mixed"])
range(all_DA_forecasts$observation[all_DA_forecasts$phen=="Stratified"])

#------------------------------------------------------------------------------#
#calculate forecast skill metrics

#forecast skill for each depth and horizon
forecast_skill_depth_horizon <-  plyr::ddply(all_DA_forecasts, c("depth","horizon", "phen", "model_id"), function(x) {
  data.frame(
        RMSE = sqrt(mean((x$mean - x$observation)^2, na.rm = TRUE)),
        MAE = mean(abs(x$mean - x$observation), na.rm = TRUE),
        pbias = 100 * (sum(x$mean - x$observation, na.rm = TRUE) / sum(x$observation, na.rm = TRUE)),
        CRPS = verification::crps(x$observation, as.matrix(x[, c(7,9)]))$CRPS,
        variance = (mean(x$sd))^2,
        sd = mean(x$sd)
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

#and last df for forecast skill aggregated across depths
forecast_avg <- plyr::ddply(all_DA_forecasts, c("model_id", "horizon"), function(x) {
  data.frame(
    RMSE = sqrt(mean((x$mean - x$observation)^2, na.rm = TRUE)),
    MAE = mean(abs(x$mean - x$observation), na.rm = TRUE),
    pbias = 100 * (sum(x$mean - x$observation, na.rm = TRUE) / sum(x$observation, na.rm = TRUE)),
    CRPS = verification::crps(x$observation, as.matrix(x[, c(7,9)]))$CRPS,
    variance = (mean(x$sd))^2
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

#order DA frequencies
forecast_avg$model_id <- factor(forecast_avg$model_id, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))


#------------------------------------------------------------------------------#
#FIGURES
cb_friendly_2 <- c("#8C510A", "#BF812D", "#C7EAE5", "#35978F") #"#DFC27D", "#DEDEDE", "#8fd5cb"


#FIGURE 5: DA frequency boxplots for 1, 5, and 9m and 1, 7, and 35-day horizons

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
#rslt_strat_9m <- c("a","c","b","c")
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
rslt_strat_1d <- c("a", "ab", "bc", "c")
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
rslt_mix_1d <- c("a","ab","bc","c")
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


#median RMSE table for depths
median_RMSE_depth <- data.frame("Depth_m" = c(rep(1,8),rep(5,8),rep(9,8)),
                                  "model_id" = rep(c("Daily","Weekly","Fortnightly","Monthly"),6),
                                  "TempDynamics" = rep(c(rep("Mixed",4),rep("Stratified",4)),3),
                                  "RMSE_C" = c(median_mix_1m_daily,median_mix_1m_weekly,median_mix_1m_fortnightly,median_mix_1m_monthly,
                                               median_strat_1m_daily,median_strat_1m_weekly,median_strat_1m_fortnightly,median_strat_1m_monthly,
                                               median_mix_5m_daily,median_mix_5m_weekly,median_mix_5m_fortnightly,median_mix_5m_monthly,
                                               median_strat_5m_daily,median_strat_5m_weekly,median_strat_5m_fortnightly,median_strat_5m_monthly,
                                               median_mix_9m_daily,median_mix_9m_weekly,median_mix_9m_fortnightly,median_mix_9m_monthly,
                                               median_strat_9m_daily,median_strat_9m_weekly,median_strat_9m_fortnightly,median_strat_9m_monthly))
#write.csv(median_RMSE_depth,file.path(lake_directory,"analysis/data/median_RMSE_depths_DA.csv"),row.names = FALSE)

#add binned RMSE for new fig 
median_RMSE_depth$RMSE_bins <- ifelse(median_RMSE_depth$RMSE_C < 1, 1, # < 1
                                        ifelse(median_RMSE_depth$RMSE_C >= 1 & median_RMSE_depth$RMSE_C < 1.5, 1.5, # 1 - 1.5
                                        ifelse(median_RMSE_depth$RMSE_C >= 1.5 & median_RMSE_depth$RMSE_C < 2, 2, # 1.5 - 2
                                               ifelse(median_RMSE_depth$RMSE_C >= 2 & median_RMSE_depth$RMSE_C < 2.5, 2.5, # 2 - 2.5
                                                      ifelse(median_RMSE_depth$RMSE_C >= 2.5, 3, 0))))) # 2.5 - 3 

#median RMSE table for horizons
median_RMSE_horizon <- data.frame("Horizon_days" = c(rep(1,8),rep(7,8),rep(35,8)),
                                  "model_id" = rep(c("Daily","Weekly","Fortnightly","Monthly"),6),
                                  "TempDynamics" = rep(c(rep("Mixed",4),rep("Stratified",4)),3),
                                  "RMSE_C" = c(median_mix_1d_daily,median_mix_1d_weekly,median_mix_1d_fortnightly,median_mix_1d_monthly,
                                               median_strat_1d_daily,median_strat_1d_weekly,median_strat_1d_fortnightly,median_strat_1d_monthly,
                                               median_mix_7d_daily,median_mix_7d_weekly,median_mix_7d_fortnightly,median_mix_7d_monthly,
                                               median_strat_7d_daily,median_strat_7d_weekly,median_strat_7d_fortnightly,median_strat_7d_monthly,
                                               median_mix_35d_daily,median_mix_35d_weekly,median_mix_35d_fortnightly,median_mix_35d_monthly,
                                               median_strat_35d_daily,median_strat_35d_weekly,median_strat_35d_fortnightly,median_strat_35d_monthly))
#write.csv(median_RMSE_horizon,file.path(lake_directory,"analysis/data/median_RMSE_horizon_DA.csv"),row.names = FALSE)

#add binned RMSE for new fig 
median_RMSE_horizon$RMSE_bins <- ifelse(median_RMSE_horizon$RMSE_C < 0.5, 0.5, # < 0.5 
                                        ifelse(median_RMSE_horizon$RMSE_C >= 0.5 & median_RMSE_horizon$RMSE_C < 1.0, 1, # 0.5 - 1 
                                               ifelse(median_RMSE_horizon$RMSE_C >= 1 & median_RMSE_horizon$RMSE_C < 1.5, 1.5, # 1 - 1.5
                                                      ifelse(median_RMSE_horizon$RMSE_C >= 1.5 & median_RMSE_horizon$RMSE_C < 2, 2, # 1.5 - 2
                                                             ifelse(median_RMSE_horizon$RMSE_C >= 2 & median_RMSE_horizon$RMSE_C < 2.5, 2.5, # 2 - 2.5
                                                                   ifelse(median_RMSE_horizon$RMSE_C >= 2.5, 3, 0)))))) # 2.5 - 3 



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

mean(median_RMSE_depth$RMSE_C[median_RMSE_depth$Depth_m==1])
mean(median_RMSE_depth$RMSE_C[median_RMSE_depth$Depth_m==5])
mean(median_RMSE_depth$RMSE_C[median_RMSE_depth$Depth_m==9])

median(median_RMSE_depth$RMSE_C[median_RMSE_depth$Depth_m==1 & median_RMSE_depth$TempDynamics=="Mixed"])
median(median_RMSE_depth$RMSE_C[median_RMSE_depth$Depth_m==5 & median_RMSE_depth$TempDynamics=="Mixed"])
median(median_RMSE_depth$RMSE_C[median_RMSE_depth$Depth_m==9 & median_RMSE_depth$TempDynamics=="Mixed"])

mean(median_RMSE_depth$RMSE_C[median_RMSE_depth$Depth_m==1 & median_RMSE_depth$TempDynamics=="Stratified"])
mean(median_RMSE_depth$RMSE_C[median_RMSE_depth$Depth_m==5 & median_RMSE_depth$TempDynamics=="Stratified"])
mean(median_RMSE_depth$RMSE_C[median_RMSE_depth$Depth_m==9 & median_RMSE_depth$TempDynamics=="Stratified"])

median(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==1 & median_RMSE_horizon$TempDynamics=="Mixed"])
median(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==7 & median_RMSE_horizon$TempDynamics=="Mixed"])
median(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==35 & median_RMSE_horizon$TempDynamics=="Mixed"])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==1 & median_RMSE_horizon$TempDynamics=="Stratified"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==7 & median_RMSE_horizon$TempDynamics=="Stratified"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==35 & median_RMSE_horizon$TempDynamics=="Stratified"])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==1 & median_RMSE_horizon$TempDynamics=="Mixed" & median_RMSE_horizon$model_id=="Daily"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==7 & median_RMSE_horizon$TempDynamics=="Mixed" & median_RMSE_horizon$model_id=="Monthly"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==35 & median_RMSE_horizon$TempDynamics=="Mixed" & median_RMSE_horizon$model_id=="Monthly"])

mean(median_RMSE_depth$RMSE_C[median_RMSE_depth$Depth_m==1 & median_RMSE_depth$TempDynamics=="Stratified" & median_RMSE_depth$model_id=="Monthly"])
mean(median_RMSE_depth$RMSE_C[median_RMSE_depth$Depth_m==5 & median_RMSE_depth$TempDynamics=="Stratified" & median_RMSE_depth$model_id=="Monthly"])
mean(median_RMSE_depth$RMSE_C[median_RMSE_depth$Depth_m==9 & median_RMSE_depth$TempDynamics=="Stratified" & median_RMSE_depth$model_id=="Fortnightly"])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==1 & median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Monthly"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==7 & median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Monthly"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==35 & median_RMSE_horizon$TempDynamics=="Stratified" & median_RMSE_horizon$model_id=="Monthly"])


#calculate 1 and 35 day RMSE for mixed vs. stratified
    mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==1 & median_RMSE_horizon$TempDynamics=="Mixed"]) 
    mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==35 & median_RMSE_horizon$TempDynamics=="Mixed"]) 

    mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==1 & median_RMSE_horizon$TempDynamics=="Stratified"]) 
    mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==35 & median_RMSE_horizon$TempDynamics=="Stratified"])
    

#mixed vs stratified rmse across depths
    mean(median_RMSE_depth$RMSE_C[median_RMSE_depth$Depth_m==1 & median_RMSE_depth$TempDynamics=="Mixed"])
    mean(median_RMSE_depth$RMSE_C[median_RMSE_depth$Depth_m==1 & median_RMSE_depth$TempDynamics=="Stratified"]) 
    
    mean(median_RMSE_depth$RMSE_C[median_RMSE_depth$Depth_m==5 & median_RMSE_depth$TempDynamics=="Mixed"])
    mean(median_RMSE_depth$RMSE_C[median_RMSE_depth$Depth_m==5 & median_RMSE_depth$TempDynamics=="Stratified"]) 
    
    mean(median_RMSE_depth$RMSE_C[median_RMSE_depth$Depth_m==9 & median_RMSE_depth$TempDynamics=="Mixed"])
    mean(median_RMSE_depth$RMSE_C[median_RMSE_depth$Depth_m==9 & median_RMSE_depth$TempDynamics=="Stratified"]) 

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#### FIGURE 5: 1,5,9m depth facets for each DA frequency for 1,7,and 35-day ahead forecasts  ####
letters <- data.frame("depth"= c(rep(1,8),rep(5,8),rep(9,8)),
                      "horizon" = rep(c(1,7,35),8),
                      "model_id" = rep(c("Daily","Weekly","Fortnightly","Monthly"),6),
                      "x" = rep(c(1,2,3,4),6),
                      "phen" = rep(c(rep("Mixed",4),rep("Stratified",4)),3),
                      "letter" = tolower(c(rslt_mix_1m, rslt_strat_1m, rslt_mix_5m, rslt_strat_5m,
                                           rslt_mix_9m,rslt_strat_9m)),
                      "max.RMSE" = c(rep(2.8,24)))

#rename depth facets
depths <- c("1m","5m","9m")
names(depths) <- c("1","5","9")

#order factor levels
kw_horizons$model_id <- factor(kw_horizons$model_id, levels = c("Daily", "Weekly", "Fortnightly", "Monthly"))

  fig_rmse_depth <- ggplot(subset(kw_horizons, horizon %in% c(1,7,35) & depth %in% c(1,5,9)),
         aes(model_id, RMSE, fill=as.factor(horizon))) +  ylab("RMSE") + xlab("")+
  geom_bar(stat="identity",position="dodge") + theme_bw() + guides(fill=guide_legend(title="Horizon")) +
  geom_hline(yintercept=2, linetype='dashed', col = 'black', linewidth=0.3)+ ylim(0,3.3) +
  scale_fill_manual(values=c("#81A665","#E0CB48","#D08151"),labels = c("1day", "7day", "35day")) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = "right",
        legend.background = element_blank(), panel.grid.minor = element_blank(), legend.box.margin=margin(-10,-1,-10,-10),
        plot.margin = unit(c(0,0.05,-0.2,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 6),legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  geom_text(data=letters,aes(x=model_id,y=0.2+max.RMSE,label=letter),hjust=0.1,vjust = -0.1, size=2.5) +
  facet_grid(depth~phen, scales="free_y",labeller = labeller(depth = depths)) 

  tag_facet2(fig_rmse_depth, fontface = 1, hjust=-0, size=3,
             tag_pool = c("a","b","c","d","e","f"))  
#ggsave(file.path(lake_directory,"analysis/figures/RMSEvsDAfreq_depth_facets.jpg"),width=3.5, height=4)

#Figure to answer Q1 - aggregated RMSE across depths and seasons
ggplot(forecast_avg ,aes(horizon, RMSE, color=as.factor(model_id))) +  
  ylab(expression("RMSE ("*~degree*C*")")) + xlab("Horizon (days)")+
  geom_line() + theme_bw() + guides(color=guide_legend(title="DA frequency")) + 
  geom_hline(yintercept = 2, linetype="dotted") + ylim(0,3.2) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = "right",
        legend.background = element_blank(),legend.direction = "vertical", panel.grid.minor = element_blank(),
        legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),legend.box.margin=margin(-10,-1,-10,-10),
        legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"), legend.title = element_text(size=6),
        axis.text.x = element_text(vjust = 0.5,size=6), axis.text.y = element_text(size=6)) +
  scale_color_manual(values=cb_friendly_2)# +geom_point()
ggsave(file.path(lake_directory,"analysis/figures/RMSEvshorizon_fig6_Q1.jpg"),width=4, height=3)

#NEW FIGURE 7!!
#now make a plot for RMSE vs all horizons
fig7 <- ggplot(subset(forecast_skill_depth_horizon, depth %in% c(1,5,9)) ,
               aes(horizon, RMSE, color=as.factor(model_id))) +  
  ylab(expression("RMSE ("*~degree*C*")")) + xlab("Horizon (days)")+
   geom_line() + theme_bw() + guides(color=guide_legend(title="DA frequency")) + 
  geom_hline(yintercept = 2, linetype="dotted") + ylim(0,3.2) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), 
        legend.position = "right", legend.background = element_blank(),
        panel.grid.minor = element_blank(), legend.key.size = unit(0.5, "lines"), 
        panel.grid.major = element_blank(),legend.box.margin=margin(-10,-1,-10,-10),
        legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"), 
        legend.title = element_text(size=6), axis.text.x = element_text(vjust = 0.5,size=6), 
        axis.text.y = element_text(size=6)) + #geom_point() +
  facet_grid(depth~phen, scales="free",labeller = labeller(depth = depths)) + 
  scale_color_manual(values=cb_friendly_2) 

tag_facet2(fig7, fontface = 1, size=3,
           tag_pool = c("a","b","c","d","e","f"))
ggsave(file.path(lake_directory,"analysis/figures/RMSEvshorizon_depth_facets_fig7.jpg"),width=3.5, height=4)


fig8 <- ggplot(subset(forecast_skill_depth_horizon, depth %in% c(1,5,9)) ,aes(horizon, variance, color=as.factor(model_id))) +  
  ylab("Variance") + xlab("Horizon (days)")+
  geom_line() + theme_bw() + guides(color=guide_legend(title="DA frequency")) + ylim(0,8.3) + 
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = "right",
        legend.background = element_blank(),legend.direction = "vertical", panel.grid.minor = element_blank(),
        legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),legend.box.margin=margin(-10,-1,-10,-10),
        legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"), legend.title = element_text(size=6),
        axis.text.x = element_text(vjust = 0.5,size=6), axis.text.y = element_text(size=6)) +
  facet_grid(depth~phen, scales="free",labeller = labeller(depth = depths)) + scale_color_manual(values=cb_friendly_2) 
tag_facet2(fig8, fontface = 1, size=3,
           tag_pool = c("a","b","c","d","e","f"))
ggsave(file.path(lake_directory,"analysis/figures/Varvshorizon_depth_facets_fig8.jpg"),width=3.5, height=4)


mean(forecast_skill_depth_horizon$RMSE)
mean(forecast_skill_depth_horizon$sd)

mean(forecast_skill_depth_horizon$RMSE[forecast_skill_depth_horizon$horizon==1])
mean(forecast_skill_depth_horizon$sd[forecast_skill_depth_horizon$horizon==1])

mean(forecast_skill_depth_horizon$RMSE[forecast_skill_depth_horizon$horizon==7])
mean(forecast_skill_depth_horizon$sd[forecast_skill_depth_horizon$horizon==7])

mean(forecast_skill_depth_horizon$RMSE[forecast_skill_depth_horizon$horizon==35])
mean(forecast_skill_depth_horizon$sd[forecast_skill_depth_horizon$horizon==35])

mean(forecast_skill_depth_horizon$RMSE[forecast_skill_depth_horizon$phen=="Mixed"])
mean(forecast_skill_depth_horizon$sd[forecast_skill_depth_horizon$phen=="Mixed"])

mean(forecast_skill_depth_horizon$RMSE[forecast_skill_depth_horizon$phen=="Stratified"])
mean(forecast_skill_depth_horizon$sd[forecast_skill_depth_horizon$phen=="Stratified"])

mean(forecast_skill_depth_horizon$RMSE[forecast_skill_depth_horizon$depth==1])
mean(forecast_skill_depth_horizon$sd[forecast_skill_depth_horizon$depth==1])

mean(forecast_skill_depth_horizon$RMSE[forecast_skill_depth_horizon$depth==5])
mean(forecast_skill_depth_horizon$sd[forecast_skill_depth_horizon$depth==5])

mean(forecast_skill_depth_horizon$RMSE[forecast_skill_depth_horizon$depth==9])
mean(forecast_skill_depth_horizon$sd[forecast_skill_depth_horizon$depth==9])

#---------------------------------------------------------------------------------#
#### SI figure: mixed and stratified RMSE tileplots aggregated across all depths  ####

#round rmse to nearest 0.5 for tile plot below
forecast_horizon_avg$RMSE_bins <- ifelse(forecast_horizon_avg$RMSE < 0.5, 0.5, # < 0.5 
                                         ifelse(forecast_horizon_avg$RMSE >= 0.5 & forecast_horizon_avg$RMSE < 1.0, 1, # 0.5 - 1 
                                                ifelse(forecast_horizon_avg$RMSE >= 1 & forecast_horizon_avg$RMSE < 1.5, 1.5, # 1 - 1.5
                                                       ifelse(forecast_horizon_avg$RMSE >= 1.5 & forecast_horizon_avg$RMSE < 2, 2, # 1.5 - 2
                                                              ifelse(forecast_horizon_avg$RMSE >= 2 & forecast_horizon_avg$RMSE < 2.5, 2.5, # 2 - 2.5
                                                                     ifelse(forecast_horizon_avg$RMSE >= 2.5, 3, 0)))))) # 2.5 - 3 

#figure for horizon vs frequency to compare forecast skill
figs8 <- ggplot(forecast_horizon_avg, aes(model_id, horizon, fill=RMSE_bins)) + 
  ylab("Horizon (days)") + theme_bw() + facet_wrap(~phen) + geom_tile(width=0.8) + xlab("") +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.spacing=unit(0, "cm"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  guides(fill=guide_legend(title="RMSE")) +  scale_fill_gradientn(labels=c(
    "< 0.5","0.5 - 1.0","1.0 - 1.5","1.5 - 2.0","2.0 - 2.5","2.5 - 3.0"),colors = hcl.colors(4, "BuPu")) 

tag_facet2(figs8, fontface = 1, hjust=0, tag_pool = c("a","b"), size=3) 
ggsave(file.path(lake_directory,"analysis/figures/HorizonvsDA_tileplot_figS8.jpg"),width=4, height=3.5)

#median mixed vs stratified period rmse across different horizons
median(forecast_horizon_avg$RMSE_bins[forecast_horizon_avg$phen=="Mixed" & forecast_horizon_avg$horizon==5 & forecast_horizon_avg$model_id=="Fortnightly"])
median(forecast_horizon_avg$RMSE_bins[forecast_horizon_avg$phen=="Stratified" & forecast_horizon_avg$horizon==2 & forecast_horizon_avg$model_id=="Weekly"])

#------------------------------------------------------------------------------------------------#
# Data assimilation figure 4
DA <- all_DA_forecasts

#change reference datetime format
DA$reference_datetime <- as.POSIXct(DA$reference_datetime)

#pull out 1 horizon for fig 4 (stratified so that all DA freqs assimilate data on same day)
all_da_sub <- DA[DA$reference_datetime=="2021-07-22" & DA$model_id=="Daily" | 
                 DA$reference_datetime=="2021-07-16" & DA$model_id=="Weekly" |
                   DA$reference_datetime=="2021-07-09" & DA$model_id=="Fortnightly" |
                   DA$reference_datetime=="2021-06-25" & DA$model_id=="Monthly",] 

#change 10m to 9m for horizon 0-2 for fortnightly example
all_da_sub$depth[all_da_sub$depth==10 & all_da_sub$model_id=="Fortnightly"] <- 9


obs <- DA[DA$reference_datetime>="2021-06-24" & DA$reference_datetime <= "2021-08-25" &
            DA$horizon==1 & DA$model_id=="Daily",] #DA for one month starting on the day that all DA frequencies assimilate data

#change date format
obs$datetime <- as.Date(obs$datetime)

#summary df to average the forecasts for each DA freq, horizon, depth, and date
all_da_sub_final <- plyr::ddply(all_da_sub, c("depth","horizon","reference_datetime", "model_id"), function(x) {
  data.frame(
    forecast_mean = mean(x$mean, na.rm = TRUE),
    forecast_sd = mean(x$sd, na.rm = TRUE),
    forecast_upper_95 = mean(x$quantile97.5, na.rm=TRUE),
    forecast_lower_95 = mean(x$quantile02.5, na.rm=TRUE),
    observed = mean(x$observation, na.rm=TRUE)
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 


#change datetime format
all_da_sub_final$reference_datetime <- as.Date(all_da_sub_final$reference_datetime)

#round depth to nearest m
all_da_sub_final$depth <- floor(all_da_sub_final$depth)
obs$depth <- floor(obs$depth)

#change DA factor order
all_da_sub_final$model_id <- factor(all_da_sub_final$model_id, levels = c("Daily", "Weekly","Fortnightly","Monthly"))

#add new column for date (reference_datetime + horizon)
all_da_sub_final$date <- all_da_sub_final$reference_datetime + all_da_sub_final$horizon

#change first 10 to 9m for 6jul-13Jul - super hacky...
obs$depth[obs$datetime=="2021-07-06" & obs$depth==10][1] <- 9
obs$depth[obs$datetime=="2021-07-07" & obs$depth==10][1] <- 9
obs$depth[obs$datetime=="2021-07-08" & obs$depth==10][1] <- 9
obs$depth[obs$datetime=="2021-07-09" & obs$depth==10][1] <- 9
obs$depth[obs$datetime=="2021-07-10" & obs$depth==10][1] <- 9
obs$depth[obs$datetime=="2021-07-11" & obs$depth==10][1] <- 9
obs$depth[obs$datetime=="2021-07-13" & obs$depth==10][1] <- 9

#change order of DA frequencies so daily is plotted on top
all_da_sub_final$model_id <- factor(all_da_sub_final$model_id, levels=rev(levels(all_da_sub_final$model_id)))

#color most recent DA date for each freq
obs$col <- ifelse(obs$datetime=="2021-07-22", "Daily",
                    ifelse(obs$datetime=="2021-07-16", "Weekly",
                           ifelse(obs$datetime=="2021-07-09", "Fortnightly",
                                  ifelse(obs$datetime=="2021-06-25", "Monthly", NA))))

obs$siz <- ifelse(is.na(obs$col), "small", "big")

#order DA freq in col column
obs$col <- factor(obs$col, levels=levels(all_da_sub_final$model_id))

ggplot(subset(all_da_sub_final, depth %in% c(1,5,9)), aes(horizon, forecast_mean, color=model_id)) + 
  geom_ribbon(data=subset(all_da_sub_final, depth %in% c(1,5,9) & model_id=="Daily"), 
              aes(x=date, y = forecast_mean, ymin = forecast_mean-forecast_sd, 
                  ymax = forecast_mean+forecast_sd), color=cb_friendly_2[1], 
              fill=cb_friendly_2[1], alpha=0.4) + 
  geom_ribbon(data=subset(all_da_sub_final, depth %in% c(1,5,9) & model_id=="Weekly"), 
              aes(x=date, y = forecast_mean, ymin = forecast_mean-forecast_sd, 
                  ymax = forecast_mean+forecast_sd), color=cb_friendly_2[2], 
              fill=cb_friendly_2[2], alpha=0.4) + 
  geom_ribbon(data=subset(all_da_sub_final, depth %in% c(1,5,9) & model_id=="Fortnightly"), 
              aes(x=date, y = forecast_mean, ymin = forecast_mean-forecast_sd, 
                  ymax = forecast_mean+forecast_sd), color=cb_friendly_2[3], 
              fill=cb_friendly_2[3], alpha=0.4) + 
  geom_ribbon(data=subset(all_da_sub_final, depth %in% c(1,5,9) & model_id=="Monthly"), 
              aes(x=date, y = forecast_mean, ymin = forecast_mean-forecast_sd, 
                  ymax = forecast_mean+forecast_sd), color=cb_friendly_2[4], 
              fill=cb_friendly_2[4], alpha=0.4) + theme_bw() + 
  scale_fill_manual(values=rev(cb_friendly_2), 
                    labels=c("Daily","Weekly","Fortnightly","Monthly")) + 
  guides(fill="none") + new_scale_fill() +
  geom_vline(xintercept=as.Date("2021-07-22")) +
  geom_point(data=subset(obs, depth %in% c(1,5,9)), 
             aes(x=datetime, y=observation, size=siz, fill=col), 
             pch=21, color="black") + 
  scale_fill_manual(values=c(rev(cb_friendly_2)), 
                    labels=c("Daily","Weekly","Fortnightly","Monthly", "")) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), 
        legend.position = "right", legend.box.margin=margin(-10,-1,-10,-10),
        legend.background = element_blank(), panel.grid.minor = element_blank(), 
        legend.key=element_rect(fill=NA), plot.margin = unit(c(0,0.05,-0.2,0), "cm"),
        legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 4.5),legend.text  = element_text(size = 4.5), 
        panel.spacing=unit(0, "cm"), axis.text.y = element_text(size=6),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6)) +
  facet_wrap(~depth, ncol = 1, scales = "free_y") +
  ylab(expression("Temperature ("*~degree*C*")")) + xlab("")  + 
  scale_color_manual(values=rev(cb_friendly_2)) + 
  scale_size_manual(values=c(0.8,0.1)) +
  geom_text(x=as.Date("2021-07-18"), y=24, label="Past", color="black", size=2) + #31.7
  geom_text(x=as.Date("2021-07-26"), y=24, label="Future", color="black", size=2) +
  guides(fill = guide_legend(title="DA frequency", 
                             override.aes = list(alpha=1,
                                                pch=c(rep(21,4),NA),
                                                fill=c(cb_friendly_2,NA)),reverse = FALSE), 
         color="none", size="none")
ggsave(file.path(lake_directory,"analysis/figures/forecast_ProofOfConcept_fig4.jpg"), width=3.5, height=4) 


#% uncertainty at 0-day horizon
mean(all_da_sub_final$forecast_upper_95[all_da_sub_final$horizon==1 & all_da_sub_final$model_id=="Daily"]) -
mean(all_da_sub_final$forecast_lower_95[all_da_sub_final$horizon==1 & all_da_sub_final$model_id=="Daily"])

mean(all_da_sub_final$forecast_upper_95[all_da_sub_final$horizon==29 & all_da_sub_final$model_id=="Monthly"]) -
  mean(all_da_sub_final$forecast_lower_95[all_da_sub_final$horizon==29 & all_da_sub_final$model_id=="Monthly"])



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

sub <- wtemp[wtemp$date > "2021-01-01" & wtemp$date < "2022-01-01", ] # Subset to target dates - need to change for your experiment
sub$Date <- as.Date(sub$date) # Line
sub$dens <- rLakeAnalyzer::water.density(sub$observation) # density for stratification

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

#rename 0 to 0.1
sub$depth[sub$depth==0] <- 0.1

#-------------------------------------------------------------------------------#
# Fig 3 - phenology plot w/ temp at different depths 
ggplot(sub) +   theme_bw() +
  geom_rect(data = phen, aes(fill = "Mixed"), xmin=-Inf ,xmax = as.Date("2021-03-11"), ymin = -Inf, ymax = Inf, inherit.aes = FALSE) + 
  geom_rect(data = phen, aes(fill = "Stratified"), xmin=as.Date("2021-03-12") ,xmax = as.Date("2021-11-07"), ymin = -Inf, ymax = Inf, inherit.aes = FALSE)+
  geom_rect(data = phen, aes(fill = "Mixed"), xmin=as.Date("2021-11-08") ,xmax = Inf, ymin = -Inf, ymax = Inf, inherit.aes = FALSE) +
  geom_line(aes(Date, as.numeric(observation), color = factor(depth)), linewidth=0.4) + ylab(expression("Temperature ("*~degree*C*")")) + xlab("") +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = "right", legend.background = element_blank(),
        legend.key = element_blank(), legend.key.height = unit(0.3,"cm"), legend.key.width = unit(0.4,"cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_fill_manual('', values = c('gray','white')) +
  scale_color_viridis(option="B", discrete="TRUE", direction=-1) +
  guides(color=guide_legend("Depth (m)"),fill= guide_legend("Period",order = 1, override.aes= list(color="black", lwd=0.1)))
ggsave(file.path(lake_directory,"analysis/figures/2021_watertemp_mixedVstratified.jpg"), width=4, height=3)

range(sub$observation[sub$datetime>= "2021-03-12" & sub$datetime<= "2021-11-07"])

#-------------------------------------------------------------------------------#
# Figure 6 - fig to compare forecast skill across horizons and depths 

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
median_RMSE_horizon$model_id <- factor(median_RMSE_horizon$model_id, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))
median_RMSE_depth$model_id <- factor(median_RMSE_depth$model_id, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))

#change horizon order
median_RMSE_horizon$Horizon_days <- factor(median_RMSE_horizon$Horizon_days, levels=c(1,7,35))
median_RMSE_depth$Depth_m <- factor(median_RMSE_depth$Depth_m, levels=c(9,5,1))

#now make a figure... 
horiz <- ggplot(median_RMSE_horizon, aes(model_id, as.factor(Horizon_days), fill=as.factor(RMSE_bins))) + 
  ylab("Horizon (days)") + xlab("") + theme_bw() + geom_tile(width=0.8) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = "right",
        legend.background = element_blank(),legend.direction = "vertical", panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0.05,0.2,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 6),legend.text  = element_text(size = 8), panel.spacing=unit(0, "cm"), legend.margin=margin(0,0,0,0),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(size=6)) +
  facet_grid(~TempDynamics, scales="free_y",labeller = labeller(Depth_m = depths)) +
  guides(fill=guide_legend(title=expression(paste("RMSE (",degree,"C)")))) + scale_fill_brewer(type="seq", 
    palette=7, direction = -1, labels = c("< 0.5", "0.5 - 1", "1 - 1.5", "1.5 - 2", "2 - 2.5", "> 2.5"))
horizon_letters <- tag_facet2(horiz, fontface = 1, hjust=0, size=3,
           tag_pool = c("a","b"))


depth <- ggplot(median_RMSE_depth, aes(model_id, as.factor(Depth_m), fill=as.factor(RMSE_bins))) + 
  ylab("Depth (m)") + xlab("") + theme_bw() + geom_tile(width=0.8) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = "none",
        legend.background = element_blank(),legend.direction = "vertical", panel.grid.minor = element_blank(),
        plot.margin = unit(c(-0.6,0.05,-0.2,0.1), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 6),legend.text  = element_text(size = 8), panel.spacing=unit(0, "cm"), legend.margin=margin(0,0,0,0),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  facet_grid(~TempDynamics, scales="free_y",labeller = labeller(Depth_m = depths)) +
   scale_fill_manual( 
    values = c("#E6550D","#FD8D3C","#FDAE6B","#FDD0A2","#FEEDDE"), labels = c(
      "0.5 - 1", "1 - 1.5", "1.5 - 2", "2 - 2.5", "> 2.5"))
depth_letters <- tag_facet2(depth, fontface = 1, hjust=0, size=3,
                              tag_pool = c("c","d"))

#combine depth and horizon plots
ggarrange(horizon_letters, depth_letters, 
          ncol = 1, nrow=2,
          common.legend = TRUE,
          legend="right") + geom_text(x=0.3, y=15, label="a")

ggsave(file.path(lake_directory,"analysis/figures/RMSE_bins_horizon_depth_facets.jpg"),width=3.5, height=4)

#-------------------------------------------------------------------------------#
#CRPS fig across mixed vs strat, depths, and horizons (fig 4 but for CRPS)

#kruskal wallis and dunn tests for 1m stratified crps across DA freq
kruskal.test(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$depth==1] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==1])
dunn_strat_1m <- dunnTest(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$depth==1] ~ 
                            kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==1])
rslt_strat_1m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_1m$res, threshold = 0.05)$Letter[c(1,4,2,3)])
rslt_strat_1m <- c("b","a","a","ab")
median_strat_1m_daily <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id=="Daily"])
median_strat_1m_weekly <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id=="Weekly"])
median_strat_1m_fortnightly <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id=="Fortnightly"])
median_strat_1m_monthly <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$depth==1 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 5m stratified crps across DA freq
kruskal.test(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$depth==5] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==5]))
dunn_strat_5m <- dunnTest(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$depth==5] ~ 
                            as.factor(kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==5]))
rslt_strat_5m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_5m$res, threshold = 0.05)$Letter[c(1,4,2,3)])
rslt_strat_5m <- c("b","a","ab","c")
median_strat_5m_daily <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$depth==5  & kw_horizons$model_id=="Daily"])
median_strat_5m_weekly <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$depth==5  & kw_horizons$model_id=="Weekly"])
median_strat_5m_fortnightly <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$depth==5  & kw_horizons$model_id=="Fortnightly"])
median_strat_5m_monthly <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$depth==5  & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 9m stratified crps across DA freq
kruskal.test(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$depth==9] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==9]))
dunn_strat_9m <- dunnTest(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$depth==9] ~ 
                            as.factor(kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==9]))
rslt_strat_9m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_9m$res, threshold = 0.05)$Letter[c(1,4,2,3)])
#rslt_strat_9m <- c("a","c","b","c")
median_strat_9m_daily <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id=="Daily"])
median_strat_9m_weekly <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id=="Weekly"])
median_strat_9m_fortnightly <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id=="Fortnightly"])
median_strat_9m_monthly <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$depth==9 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 1m mixed crps across DA freq
kruskal.test(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$depth==1] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==1]))
dunn_mix_1m <- dunnTest(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$depth==1] ~ 
                          as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==1]))
rslt_mix_1m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_1m$res, threshold = 0.05)$Letter[c(1,4,2,3)])
rslt_mix_1m <- c("b","a","ab","ab")
median_mix_1m_daily <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id=="Daily"])
median_mix_1m_weekly <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id=="Weekly"])
median_mix_1m_fortnightly <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id=="Fortnightly"])
median_mix_1m_monthly <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$depth==1 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 5m mixed crps across DA freq
kruskal.test(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$depth==5] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==5]))
dunn_mix_5m <- dunnTest(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$depth==5] ~ 
                          as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==5]))
rslt_mix_5m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_5m$res, threshold = 0.05)$Letter[c(1,4,2,3)])
rslt_mix_5m <- c("b","a","ab","ab")
median_mix_5m_daily <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id=="Daily"])
median_mix_5m_weekly <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id=="Weekly"])
median_mix_5m_fortnightly <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id=="Fortnightly"])
median_mix_5m_monthly <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$depth==5 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 9m mixed crps across DA freq
kruskal.test(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$depth==9] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==9]))
dunn_mix_9m <- dunnTest(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$depth==9] ~ 
                          as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$depth==9]))
rslt_mix_9m=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_9m$res, threshold = 0.05)$Letter[c(1,4,2,3)])
rslt_mix_9m <- c("b","a","ab","ab")
median_mix_9m_daily <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id=="Daily"])
median_mix_9m_weekly <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id=="Weekly"])
median_mix_9m_fortnightly <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id=="Fortnightly"])
median_mix_9m_monthly <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$depth==9 & kw_horizons$model_id=="Monthly"])

#Now aggregate across depths
#kruskal wallis and dunn tests for 1d stratified crps across DA freq
kruskal.test(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1])
dunn_strat_1d <- dunnTest(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1] ~ 
                            kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$horizon==1])
rslt_strat_1d=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_1d$res, threshold = 0.05)$Letter[c(1,4,2,3)])
rslt_strat_1d <- c("a","ab","bc","c")
median_strat_1d_daily <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==1 & kw_horizons$model_id=="Daily"])
median_strat_1d_weekly <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==1 & kw_horizons$model_id=="Weekly"])
median_strat_1d_fortnightly <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==1 & kw_horizons$model_id=="Fortnightly"])
median_strat_1d_monthly <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==1 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 7d stratified crps across DA freq
kruskal.test(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7]))
dunn_strat_7d <- dunnTest(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7] ~ 
                            as.factor(kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$horizon==7]))
rslt_strat_7d=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_7d$res, threshold = 0.05)$Letter[c(1,4,2,3)])
median_strat_7d_daily <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==7 & kw_horizons$model_id=="Daily"])
median_strat_7d_weekly <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==7 & kw_horizons$model_id=="Weekly"])
median_strat_7d_fortnightly <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==7 & kw_horizons$model_id=="Fortnightly"])
median_strat_7d_monthly <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==7 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 35d stratified crps across DA freq
kruskal.test(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35]))
dunn_strat_35d <- dunnTest(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35] ~ 
                             as.factor(kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$horizon==35]))
rslt_strat_35d=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_35d$res, threshold = 0.05)$Letter[c(1,4,2,3)])
median_strat_35d_daily <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==35 & kw_horizons$model_id=="Daily"])
median_strat_35d_weekly <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==35 & kw_horizons$model_id=="Weekly"])
median_strat_35d_fortnightly <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==35 & kw_horizons$model_id=="Fortnightly"])
median_strat_35d_monthly <- median(kw_horizons$CRPS[kw_horizons$phen=="Stratified" & kw_horizons$horizon ==35 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 1d mixed crps across DA freq
kruskal.test(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1]))
dunn_mix_1d <- dunnTest(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1] ~ 
                          as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$horizon==1]))
rslt_mix_1d=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_1d$res, threshold = 0.05)$Letter[c(1,4,2,3)])
rslt_mix_1d <- c("a","ab","bc","c")
median_mix_1d_daily <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==1 & kw_horizons$model_id=="Daily"])
median_mix_1d_weekly <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==1 & kw_horizons$model_id=="Weekly"])
median_mix_1d_fortnightly <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==1 & kw_horizons$model_id=="Fortnightly"])
median_mix_1d_monthly <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==1 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 7d mixed crps across DA freq
kruskal.test(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7]))
dunn_mix_7d <- dunnTest(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7] ~ 
                          as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$horizon==7]))
rslt_mix_7d=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_7d$res, threshold = 0.05)$Letter[c(1,4,2,3)])
rslt_mix_7d <- c("ab","a","bc","c")
median_mix_7d_daily <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==7 & kw_horizons$model_id=="Daily"])
median_mix_7d_weekly <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==7 & kw_horizons$model_id=="Weekly"])
median_mix_7d_fortnightly <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==7 & kw_horizons$model_id=="Fortnightly"])
median_mix_7d_monthly <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==7 & kw_horizons$model_id=="Monthly"])

#kruskal wallis and dunn tests for 35d mixed crps across DA freq
kruskal.test(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35] ~ 
               as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35]))
dunn_mix_35d <- dunnTest(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35] ~ 
                           as.factor(kw_horizons$model_id[kw_horizons$phen=="Mixed" & kw_horizons$horizon==35]))
rslt_mix_35d=toupper(cldList(P.adj ~ Comparison, data=dunn_mix_35d$res, threshold = 0.05)$Letter[c(1,4,2,3)])
rslt_mix_35d <- c("c","a","bc","ab")
median_mix_35d_daily <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==35 & kw_horizons$model_id=="Daily"])
median_mix_35d_weekly <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==35 & kw_horizons$model_id=="Weekly"])
median_mix_35d_fortnightly <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==35 & kw_horizons$model_id=="Fortnightly"])
median_mix_35d_monthly <- median(kw_horizons$CRPS[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==35 & kw_horizons$model_id=="Monthly"])

#NOW PLOT 
letters <- data.frame("depth"= c(rep(1,8),rep(5,8),rep(9,8)),
                      "horizon" = rep(c(1,7,35),8),
                      "model_id" = rep(c("Daily","Weekly","Fortnightly","Monthly"),6),
                      "x" = rep(c(1,2,3,4),6),
                      "phen" = rep(c(rep("Mixed",4),rep("Stratified",4)),3),
                      "letter" = tolower(c(rslt_mix_1m, rslt_strat_1m, rslt_mix_5m, rslt_strat_5m,
                                           rslt_mix_9m,rslt_strat_9m)),
                      "max.RMSE" = c(rep(2.6,8),rep(2.5,8),rep(2,8)))

#rename depth facets
depths <- c("1m","5m","9m")
names(depths) <- c("1","5","9")

#order factor levels
kw_horizons$model_id <- factor(kw_horizons$model_id, levels = c("Daily", "Weekly", "Fortnightly", "Monthly"))

figs1 <- ggplot(subset(forecast_skill_depth_horizon, depth %in% c(1,5,9)) ,
               aes(horizon, CRPS, color=as.factor(model_id))) +  
  ylab("CRPS") + xlab("Horizon (days)")+
  geom_line() + theme_bw() + guides(color=guide_legend(title="DA frequency")) + 
  ylim(0,2) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), 
        legend.position = "right", legend.background = element_blank(),
        panel.grid.minor = element_blank(), legend.key.size = unit(0.5, "lines"), 
        panel.grid.major = element_blank(),legend.box.margin=margin(-10,-1,-10,-10),
        legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"), 
        legend.title = element_text(size=6), axis.text.x = element_text(vjust = 0.5,size=6), 
        axis.text.y = element_text(size=6)) + #geom_point() +
  facet_grid(depth~phen, scales="free",labeller = labeller(depth = depths)) + 
  scale_color_manual(values=cb_friendly_2) 

tag_facet2(figs1, fontface = 1, size=3,
           tag_pool = c("a","b","c","d","e","f"))

ggsave(file.path(lake_directory,"analysis/figures/CRPSvsDAfreq_depth_facets_figs1.jpg"),width=3.5, height=4)

#------------------------------------------------------------------------------#
#figure to visualize different DA frequencies over 1 year period
daily <- seq.Date(as.Date("2020-11-27"), as.Date("2021-12-31"), by = 1) 
weekly <- daily[seq(1, length(daily), 7)] 
fortnightly <-  daily[seq(1, length(daily), 14)]
monthly <- daily[seq(1, length(daily), 30)]

daily <- data.frame(daily) %>% rename(day=daily)
weekly <- data.frame(weekly) %>% rename(day=weekly)
fortnightly <- data.frame(fortnightly) %>% rename(day=fortnightly)
monthly <- data.frame(monthly) %>% rename(day=monthly)

weekly$weekly <- seq.Date(as.Date("2020-11-27"), as.Date("2021-12-31"), by = 7) 
fortnightly$fortnightly <- seq.Date(as.Date("2020-11-27"), as.Date("2021-12-31"), by = 14) 
monthly$monthly <- seq.Date(as.Date("2020-11-27"), as.Date("2021-12-31"), by = 30) 

d <- data.frame("Daily"= seq.Date(as.Date("2020-11-27"), as.Date("2021-12-31"), by = 1) )
w<- pad(weekly, interval="day", by="day")
f <- pad(fortnightly, interval="day", by="day", end_val = as.Date("2021-12-31"))
m <- pad(monthly, interval="day", by="day", end_val = as.Date("2021-12-31"))

all_dates <- data.frame("day" = daily,"Daily"= d, "Weekly"= w$weekly,
                        "Fortnightly" = f$fortnightly, "Monthly" = m$monthly)

#wide to long
dates_long <-   all_dates %>%
  gather(model_id,value, Daily:Monthly) 

#order DA frequencies
dates_long$model_id <- factor(dates_long$model_id, levels=c("Daily", "Weekly", "Fortnightly", "Monthly")) 

#plot DA frequencies 
ggplot(dates_long, aes(value,model_id, color=model_id)) + geom_count(show.legend=FALSE, size=4, pch=124) +
  scale_color_manual(values=cb_friendly_2) + ylab("DA frequency") + xlab("") +
  theme_bw() +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(vjust = 0.5, hjust=1,size=6), 
        axis.text.y = element_text(size=6)) 
ggsave(file.path(lake_directory,"analysis/figures/DA_freq_2021_SI.jpg"), width=4, height=3.5) 

#-------------------------------------------------------------------------------#
# NEW FIGURE - RMSE vs st dev
fig_test <- ggplot(subset(forecast_skill_depth_horizon, depth %in% c(1,5,9) & horizon %in% c(1,7,35)))+
                   geom_line(aes(sd,RMSE, color=as.factor(model_id))) + geom_point(aes(sd, RMSE, color=as.factor(model_id), shape=as.factor(horizon)), size=2.5) +  ylim(0,3.3) +
  ylab(expression("RMSE ("*~degree*C*")")) + xlab(expression("Standard deviation ("*~degree*C*")"))+
   theme_bw() + guides(color=guide_legend(title="DA frequency"), shape=guide_legend(title="Horizon (days)")) + 
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = "right",
        legend.background = element_blank(),legend.direction = "vertical", panel.grid.minor = element_blank(),
        legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),legend.box.margin=margin(-10,-1,-10,-10),
        legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"), legend.title = element_text(size=6),
        axis.text.x = element_text(vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  facet_grid(depth~phen, scales="free",labeller = labeller(depth = depths)) + scale_color_manual(values=cb_friendly_2) 

tag_facet2(fig_test, fontface = 1, size=3,
           tag_pool = c("a","b","c","d","e","f"))
ggsave(file.path(lake_directory,"analysis/figures/RMSEvsStdev_depth_facets.jpg"),width=3.5, height=4)


