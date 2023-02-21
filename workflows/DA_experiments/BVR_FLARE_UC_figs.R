#UC analysis figs
#29 Dec 2022

#load libraries
pacman::p_load(ggplot2,tidyverse,FSA,rcompanion,ggpubr, egg)

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
score_dir <- arrow::SubTreeFileSystem$create(file.path(lake_directory,"scores/UC"))
all_DA_forecasts <- arrow::open_dataset(score_dir) |> collect() |>   
  filter(!is.na(observation), variable == "temperature",horizon > 1) #need to start on day 2 because this is actually day 1


#need to round horizon becuase they are in decimal form for some reason...
all_DA_forecasts$horizon <- ceiling(all_DA_forecasts$horizon)

#round depths up to nearest m 
all_DA_forecasts$depth <- ceiling(all_DA_forecasts$depth)

#fix horizon issue because flare calls day 0 day 1 (so horizons go out to 36 days)
all_DA_forecasts$horizon <- all_DA_forecasts$horizon - 1

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

#------------------------------------------------------------------------------#
#calculate forecast skill metrics

#forecast skill for each depth and horizon
forecast_skill_depth_horizon <-  plyr::ddply(all_DA_forecasts, c("depth","horizon","phen", "model_id"), function(x) { #phen instead of datetime?
  data.frame(
    RMSE = sqrt(mean((x$mean - x$observation)^2, na.rm = TRUE)),
    MAE = mean(abs(x$mean - x$observation), na.rm = TRUE),
    pbias = 100 * (sum(x$mean - x$observation, na.rm = TRUE) / sum(x$observation, na.rm = TRUE)),
    CRPS = verification::crps(x$observation, as.matrix(x[, c(7,9)]))$CRPS,
    variance = (mean(x$sd))^2,
    sd = mean(x$sd)
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 


#quickly looking at variance among dates in forecast period
#add month column
all_DA_forecasts$month <- format(as.Date(all_DA_forecasts$reference_datetime),"%b")

variance_2021_months <-  plyr::ddply(all_DA_forecasts, c("month","horizon", "model_id","depth"), function(x) { #phen instead of datetime?
  data.frame(
    mean = mean(x$mean),
    sd = mean(x$sd),
    variance = (mean(x$sd))^2
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 


ggplot(subset(variance_2021_months, model_id=="Daily" & depth==1), aes(horizon, variance)) + geom_line() +
  theme_bw() + facet_wrap(~month)


#averaging across all months
variance_2021_days_year <-  plyr::ddply(all_DA_forecasts, c("horizon", "model_id"), function(x) { #phen instead of datetime?
  data.frame(
    mean = mean(x$mean),
    sd = mean(x$sd),
    variance = (mean(x$sd))^2
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

ggplot(subset(variance_2021_days_year, model_id=="Daily"), aes(horizon, variance)) +
  geom_line() + theme_bw() 

variance_2021_days <-  plyr::ddply(all_DA_forecasts, c("datetime","horizon", "model_id","depth"), function(x) { #phen instead of datetime?
  data.frame(
    mean = mean(x$mean),
    sd = mean(x$sd),
    variance = (mean(x$sd))^2
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

#only select april

variance_2021_days$datetime <- as.Date(variance_2021_days$datetime)
variance_2021_april <- variance_2021_days[variance_2021_days$datetime %in% 
                                            seq(as.Date("2021-04-01"),as.Date("2021-04-30"),by="day"),]


#look at variance for every day in april to see which noaa forecasts may be problematic 
ggplot(subset(variance_2021_april, model_id=="Daily" & depth==1),
       aes(horizon, variance)) + geom_vline(xintercept=16, linetype=2) +
  geom_line() + theme_bw()+facet_wrap(~datetime)


#df with averaged forecast skill for all days (group by horizon, DA, and phen)
forecast_horizon_avg <- plyr::ddply(all_DA_forecasts, c("horizon", "model_id", "phen"), function(x) {
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
forecast_horizon_avg$model_id <- factor(forecast_horizon_avg$model_id, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))

#df with averaged forecast skill for all days (group by depth, DA, and phen)
forecast_depth_avg <- plyr::ddply(all_DA_forecasts, c("depth", "model_id", "phen"), function(x) {
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
forecast_depth_avg$model_id <- factor(forecast_depth_avg$model_id, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))


#------------------------------------------------------------------------------#
#FIGURES
cb_friendly_2 <- c("#8C510A", "#BF812D", "#C7EAE5", "#35978F") #"#DFC27D", "#DEDEDE",

#FIGURE 5: DA frequency boxplots for 1, 5, and 9m and 1, 7, and 35-day horizons

#creating smaller dataset for kw test w/ 1,5,9m and 1,7,35 days
kw_horizons <- forecast_skill_depth_horizon #[forecast_skill_depth_horizon$depth %in% c(1,5,9) & forecast_skill_depth_horizon$horizon %in% c(1,7,35),]

#kruskal wallis and dunn tests for 1m stratified rmse across DA freq
kruskal.test(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1] ~ 
               kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==1])
dunn_strat_1m <- dunnTest(kw_horizons$RMSE[kw_horizons$phen=="Stratified" & kw_horizons$depth==1] ~ 
                            kw_horizons$model_id[kw_horizons$phen=="Stratified" & kw_horizons$depth==1])
rslt_strat_1m=toupper(cldList(P.adj ~ Comparison, data=dunn_strat_1m$res, threshold = 0.05)$Letter[c(1,4,2,3)])
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
rslt_strat_5m <- c("b","a","b","c")
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
rslt_strat_9m <- c("a","b","c","d")
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
rslt_mix_1m <- c("c","a","ab","bc")
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
rslt_mix_5m <- c("c","a","b","bc")
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
rslt_mix_9m <- c("c","a","b","bc")
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
rslt_strat_1d <- c("a","ab","bc","c")
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
rslt_mix_7d <- c("bc","a","ab","c")
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
rslt_mix_35d <- c("c","a","b","b")
median_mix_35d_daily <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==35 & kw_horizons$model_id=="Daily"])
median_mix_35d_weekly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==35 & kw_horizons$model_id=="Weekly"])
median_mix_35d_fortnightly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==35 & kw_horizons$model_id=="Fortnightly"])
median_mix_35d_monthly <- median(kw_horizons$RMSE[kw_horizons$phen=="Mixed" & kw_horizons$horizon ==35 & kw_horizons$model_id=="Monthly"])

#create table to export mixed and stratified p-vals across depths (aggregated over horizon)
pvals_horizon_aggregated <- data.frame("Depth_m"= c(rep(1,6),rep(5,6),rep(9,6)), 
                                       "Comparison" = c(dunn_mix_1m$res$Comparison,dunn_mix_5m$res$Comparison,dunn_mix_9m$res$Comparison),
                                       #"Z" = c(dunn_mix_1m$res$Z,dunn_mix_5m$res$Z,dunn_mix_9m$res$Z),
                                       "Mixed_pvalue" = c(dunn_mix_1m$res$P.adj,dunn_mix_5m$res$P.adj,dunn_mix_9m$res$P.adj),
                                       "Stratified_pvalue" = c(dunn_strat_1m$res$P.adj,dunn_strat_5m$res$P.adj,dunn_strat_9m$res$P.adj))

#write.csv(pvals_horizon_aggregated,file.path(lake_directory,"analysis/data/UC_pvals_horizon_aggregatred.csv"),row.names = FALSE)


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
#write.csv(median_RMSE_depth,file.path(lake_directory,"analysis/data/UC_median_RMSE_depths_DA.csv"),row.names = FALSE)

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
#write.csv(median_RMSE_horizon,file.path(lake_directory,"analysis/data/UC_median_RMSE_horizon_DA.csv"),row.names = FALSE)

#add binned RMSE for new fig 
median_RMSE_horizon$RMSE_bins <- ifelse(median_RMSE_horizon$RMSE_C < 0.5, 0.5, # < 0.5 
                                        ifelse(median_RMSE_horizon$RMSE_C >= 0.5 & median_RMSE_horizon$RMSE_C < 1.0, 1, # 0.5 - 1 
                                               ifelse(median_RMSE_horizon$RMSE_C >= 1 & median_RMSE_horizon$RMSE_C < 1.5, 1.5, # 1 - 1.5
                                                      ifelse(median_RMSE_horizon$RMSE_C >= 1.5 & median_RMSE_horizon$RMSE_C < 2, 2, # 1.5 - 2
                                                             ifelse(median_RMSE_horizon$RMSE_C >= 2 & median_RMSE_horizon$RMSE_C < 2.5, 2.5, # 2 - 2.5
                                                                    ifelse(median_RMSE_horizon$RMSE_C >= 2.5, 3, 0)))))) # 2.5 - 3 


mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$TempDynamics=="Mixed"])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$TempDynamics=="Stratified"])

mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==1])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==7])
mean(median_RMSE_horizon$RMSE_C[median_RMSE_horizon$Horizon_days==35])

mean(median_RMSE_depth$RMSE_C[median_RMSE_depth$Depth_m==1])
mean(median_RMSE_depth$RMSE_C[median_RMSE_depth$Depth_m==5])
mean(median_RMSE_depth$RMSE_C[median_RMSE_depth$Depth_m==9])

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

ggplot(subset(kw_horizons, horizon %in% c(1,7,35) & depth %in% c(1,5,9)), 
       aes(model_id, RMSE, fill=as.factor(horizon))) +  ylab("RMSE") + xlab("")+
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

depth <- ggplot(median_RMSE_depth, aes(model_id, as.factor(Depth_m), fill=as.factor(RMSE_bins))) + 
  ylab("Depth (m)") + xlab("") + theme_bw() + geom_tile(width=0.8) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = "right",
        legend.background = element_blank(),legend.direction = "vertical", panel.grid.minor = element_blank(),
        plot.margin = unit(c(-0.6,0.05,-0.2,0.1), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 6),legend.text  = element_text(size = 8), panel.spacing=unit(0, "cm"), legend.margin=margin(0,0,0,0),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  facet_grid(~TempDynamics, scales="free_y",labeller = labeller(Depth_m = depths)) +
  guides(fill=guide_legend(title=expression(paste("RMSE (",degree,"C)")))) + scale_fill_manual( 
    values = c("#E6550D","#FD8D3C","#FDAE6B","#FDD0A2","#FEEDDE"), labels = c(
      "0.5 - 1", "1 - 1.5", "1.5 - 2", "2 - 2.5", "> 2.5"))

#combine depth and horizon plots
ggarrange(horiz, depth, 
          ncol = 1, nrow = 2,
          common.legend = TRUE,
          legend="right")


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

last(params$mean[params$variable=="lw_factor" & params$model_id=="Daily"])

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
  filter(!is.na(observation), variable == "temperature",horizon > 1)

#round depths up to nearest m 
all_DA_forecasts_yesIC$depth <- ceiling(all_DA_forecasts_yesIC$depth)

#fix horizon issue because flare calls day 0 day 1 (so horizons go out to 36 days)
all_DA_forecasts_yesIC$horizon <- all_DA_forecasts_yesIC$horizon - 1

#add a group number so that I can average horizons later on
all_DA_forecasts_yesIC <- all_DA_forecasts_yesIC %>% 
  mutate(group = case_when(all_DA_forecasts_yesIC$horizon <= 5 ~ "1-5",
                           all_DA_forecasts_yesIC$horizon <=10 & all_DA_forecasts_yesIC$horizon > 5 ~ "6-10",
                           all_DA_forecasts_yesIC$horizon <=15 & all_DA_forecasts_yesIC$horizon > 10 ~ "11-15",
                           all_DA_forecasts_yesIC$horizon <=20 & all_DA_forecasts_yesIC$horizon > 15 ~ "16-20",
                           all_DA_forecasts_yesIC$horizon <=25 & all_DA_forecasts_yesIC$horizon > 20 ~ "21-25",
                           all_DA_forecasts_yesIC$horizon <=30 & all_DA_forecasts_yesIC$horizon > 25 ~ "26-30",
                           all_DA_forecasts_yesIC$horizon <=35 & all_DA_forecasts_yesIC$horizon > 30 ~ "31-35"))

strat_date<- "2021-11-07"

#add stratified vs mixed col
all_DA_forecasts_yesIC$phen <- ifelse(all_DA_forecasts_yesIC$datetime <= as.POSIXct(strat_date) & 
                                        all_DA_forecasts_yesIC$datetime >="2021-03-12","Stratified", "Mixed")

#remove n=6 days with ice-cover 
all_DA_forecasts_yesIC <- all_DA_forecasts_yesIC[!(as.Date(all_DA_forecasts_yesIC$datetime) %in% c(as.Date("2021-01-10"), as.Date("2021-01-11"),as.Date("2021-01-30"),
                                                                                                   as.Date("2021-02-13"),as.Date("2021-02-14"),as.Date("2021-02-15"))),]
#drop 11m completely because some rows were NA when water level was low
all_DA_forecasts_yesIC <- all_DA_forecasts_yesIC[!(all_DA_forecasts_yesIC$depth==11),]

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
    variance = (mean(x$sd))^2,
    sd = mean(x$sd)
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
    CRPS = verification::crps(x$observation, as.matrix(x[, c(7,9)]))$CRPS,
    variance = (mean(x$sd))^2,
    sd = mean(x$sd)
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
    variance = (mean(x$sd))^2,
    sd = mean(x$sd)
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

IC <- ggplot(subset(UC_depth, depth %in% c(1,5,9)), aes(model_id, RMSE, fill=as.factor(IC))) +  ylab("RMSE") + xlab("")+
  geom_bar(stat="identity",position="dodge") + theme_bw() + guides(fill=guide_legend(title="IC")) + geom_hline(yintercept=2,linetype="dashed") +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = "right",
        legend.background = element_blank(), panel.grid.minor = element_blank(), legend.box.margin=margin(-10,-1,-10,-10),
        plot.margin = unit(c(0,0.05,-0.2,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 6),legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  facet_grid(depth~phen, scales="free",labeller = labeller(depth = depths)) + scale_fill_manual(values=c("#81A665","#E0CB48")) 

tag_facet2(IC, fontface = 1, hjust=0, size=3, tag_pool = c("a","b","c","d","e","f"))
ggsave(file.path(lake_directory,"analysis/figures/UC_RMSEvsDAfreq_depth_facets_IC_allhorizons.jpg"),width=3.5, height=4)

IC_var <- ggplot(subset(UC, depth %in% c(1,5,9) & horizon==1), aes(model_id, sd, fill=as.factor(IC))) +  ylab("sd") + xlab("")+
  geom_bar(stat="identity",position="dodge") + theme_bw() + guides(fill=guide_legend(title="")) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = "right",
        legend.background = element_blank(), panel.grid.minor = element_blank(), legend.box.margin=margin(-10,-1,-10,-10),
        plot.margin = unit(c(0,0.05,-0.2,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 6),legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  facet_grid(depth~phen, scales="free",labeller = labeller(depth = depths)) + scale_fill_manual(values=c("#81A665","#E0CB48")) 

tag_facet2(IC_var, fontface = 1, hjust=0, size=3, tag_pool = c("a","b","c","d","e","f"))
ggsave(file.path(lake_directory,"analysis/figures/UC_SDevvsDAfreq_depth_facets_IC_1day.jpg"),width=3.5, height=4)


fig8 <- ggplot(subset(UC_depth, depth %in% c(1,5,9)) ,aes(model_id, variance, fill=as.factor(IC))) +  ylab("variance") + xlab("")+
  geom_bar(stat="identity",position="dodge") + theme_bw() + guides(fill=guide_legend(title="IC")) + 
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = "right",
        legend.background = element_blank(), panel.grid.minor = element_blank(), legend.box.margin=margin(-10,-1,-10,-10),
        plot.margin = unit(c(0,0.05,-0.2,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 6),legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  facet_grid(depth~phen, scales="free",labeller = labeller(depth = depths)) + scale_fill_manual(values=c("#81A665","#E0CB48")) 

tag_facet2(fig8, fontface = 1, size=3, tag_pool = c("a","b","c","d","e","f"))  
ggsave(file.path(lake_directory,"analysis/figures/UC_variancevsDAfreq_depth_facets_IC.jpg"),width=3.5, height=4)

mean(UC_depth$variance[UC_depth$phen=="Mixed" & UC_depth$IC=="yes IC"])
mean(UC_depth$variance[UC_depth$phen=="Mixed" & UC_depth$IC=="no IC"])

sd(UC_depth$variance[UC_depth$phen=="Mixed" & UC_depth$IC=="yes IC"])
sd(UC_depth$variance[UC_depth$phen=="Mixed" & UC_depth$IC=="no IC"])

mean(UC_depth$variance[UC_depth$phen=="Stratified" & UC_depth$IC=="yes IC" & UC_depth$depth==1])
mean(UC_depth$variance[UC_depth$phen=="Stratified" & UC_depth$IC=="no IC"& UC_depth$depth==1])

UC_depth$variance[UC_depth$phen=="Stratified" & UC_depth$IC=="yes IC" & UC_depth$depth==5 & UC_depth$model_id=="Monthly"]
UC_depth$variance[UC_depth$phen=="Stratified" & UC_depth$IC=="no IC"& UC_depth$depth==5 & UC_depth$model_id=="Monthly"]

UC_depth$variance[UC_depth$phen=="Stratified" & UC_depth$IC=="yes IC" & UC_depth$depth==9 & UC_depth$model_id=="Monthly"]
UC_depth$variance[UC_depth$phen=="Stratified" & UC_depth$IC=="no IC"& UC_depth$depth==9 & UC_depth$model_id=="Monthly"]

#-------------------------------------------------------------------------------#
#NEW FIGURE - proportion of IC uncertainty vs horizon

#new dataframe to calculate the proportion of IC uncertainty across depths, horizons, model_id, and phen
UC_prop <- data.frame("depth" = forecast_skill_depth_horizon_yesIC$depth, 
                      "horizon" = forecast_skill_depth_horizon_yesIC$horizon, 
                      "phen" = forecast_skill_depth_horizon_yesIC$phen, 
                      "model_id" = forecast_skill_depth_horizon_yesIC$model_id,
                      "var" = forecast_skill_depth_horizon_yesIC$variance, 
                      "sd" = forecast_skill_depth_horizon_yesIC$sd)

UC_prop$prop_var <- ((forecast_skill_depth_horizon_yesIC$variance - 
                        forecast_skill_depth_horizon$variance) / 
                       forecast_skill_depth_horizon_yesIC$variance)

UC_prop$prop_sd <- ((forecast_skill_depth_horizon_yesIC$sd - 
                       forecast_skill_depth_horizon$sd) / 
                      forecast_skill_depth_horizon_yesIC$sd)


fig9 <- ggplot(subset(UC_prop, depth %in% c(1,5,9)) ,
               aes(horizon, prop_sd, color=as.factor(model_id))) +  ylab("Proportion of IC uncertainty") + 
  geom_line() + theme_bw() + guides(color=guide_legend(title="DA frequency")) + xlab("Horizon (Days)")+
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = "right",
        legend.background = element_blank(),legend.direction = "vertical", panel.grid.minor = element_blank(),
        legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),legend.box.margin=margin(-10,-1,-10,-10),
        legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"), legend.title = element_text(size=6),
        axis.text.x = element_text(vjust = 0.5,size=6), axis.text.y = element_text(size=6)) + ylim(-0.06, 0.75) +
  facet_grid(depth~phen, scales="free",labeller = labeller(depth = depths)) + scale_color_manual(values=cb_friendly_2) 

tag_facet2(fig9, fontface = 1, size=3, tag_pool = c("a","b","c","d","e","f"))  
ggsave(file.path(lake_directory,"analysis/figures/UC_prop_ICdvshorizon_depth_facets_sd.jpg"),width=3.5, height=4)

mean(UC_prop$prop_sd[UC_prop$horizon==1 & UC_prop$model_id=="Fortnightly"])

mean(UC_prop$prop_sd[UC_prop$phen=="Mixed"])
mean(UC_prop$prop_sd[UC_prop$phen=="Stratified" & UC_prop$depth==9])

