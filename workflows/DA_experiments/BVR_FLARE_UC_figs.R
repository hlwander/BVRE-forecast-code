library(tidyverse)
library(lubridate)
#UC analysis figs
#29 Dec 2022

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
score_dir <- arrow::SubTreeFileSystem$create(file.path(lake_directory,"scores/IC_off"))
all_DA_forecasts_noIC <- arrow::open_dataset(score_dir) |>    
  filter(!is.na(observation), variable == "temperature",horizon > 0.3) |>
  collect() 


#need to round horizon because they are in decimal form for all Jan 1 ref date
all_DA_forecasts_noIC$horizon <- floor(all_DA_forecasts_noIC$horizon)

#round depths up to nearest m 
#all_DA_forecasts$depth <- ceiling(all_DA_forecasts$depth)

#add a group number so that I can average horizons later on
all_DA_forecasts_noIC <- all_DA_forecasts_noIC %>% 
  mutate(group = case_when(all_DA_forecasts_noIC$horizon <= 5 ~ "1-5",
                           all_DA_forecasts_noIC$horizon <=10 & all_DA_forecasts_noIC$horizon > 5 ~ "6-10",
                           all_DA_forecasts_noIC$horizon <=15 & all_DA_forecasts_noIC$horizon > 10 ~ "11-15",
                           all_DA_forecasts_noIC$horizon <=20 & all_DA_forecasts_noIC$horizon > 15 ~ "16-20",
                           all_DA_forecasts_noIC$horizon <=25 & all_DA_forecasts_noIC$horizon > 20 ~ "21-25",
                           all_DA_forecasts_noIC$horizon <=30 & all_DA_forecasts_noIC$horizon > 25 ~ "26-30",
                           all_DA_forecasts_noIC$horizon <=35 & all_DA_forecasts_noIC$horizon > 30 ~ "31-35"))

strat_date<- "2021-11-07"

#add stratified vs mixed col
all_DA_forecasts_noIC$phen <- ifelse(all_DA_forecasts_noIC$datetime <= as.POSIXct(strat_date) & 
                                       all_DA_forecasts_noIC$datetime >="2021-03-12","Stratified", "Mixed")

#remove n=6 days with ice-cover 
all_DA_forecasts_noIC <- all_DA_forecasts_noIC[!(as.Date(all_DA_forecasts_noIC$datetime) %in% c(as.Date("2021-01-10"), as.Date("2021-01-11"),as.Date("2021-01-30"),
                                                                                 as.Date("2021-02-13"),as.Date("2021-02-14"),as.Date("2021-02-15"))),]

#drop 11m completely because some rows were NA when water level was low
all_DA_forecasts_noIC <- all_DA_forecasts_noIC[!(all_DA_forecasts_noIC$depth==11),]

#change model_id to be all uppercase
all_DA_forecasts_noIC$model_id <- str_to_title(all_DA_forecasts_noIC$model_id)

#only keep 2021 data
all_DA_forecasts_noIC <- all_DA_forecasts_noIC[all_DA_forecasts_noIC$datetime<="2021-12-31",]

#------------------------------------------------------------------------------#
#calculate forecast skill metrics

#add new column for binning different layers of the water column
all_DA_forecasts_noIC$layer <- ifelse(
  all_DA_forecasts_noIC$depth <= 2.5, "surface",
  ifelse(all_DA_forecasts_noIC$depth >= 8.5, "bottom", "middle"))

#forecast skill for each depth and horizon
forecast_skill_layer_horizon_noIC <-  plyr::ddply(all_DA_forecasts_noIC, c("layer","horizon","phen", "model_id"), function(x) { #phen instead of datetime?
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
forecast_skill_layer_horizon_noIC$model_id <- 
  factor(forecast_skill_layer_horizon_noIC$model_id, levels=c(
    "Daily", "Weekly", "Fortnightly", "Monthly"))


#df with averaged forecast skill for all days (group by horizon, DA, and phen)
forecast_horizon_avg_noIC <- plyr::ddply(all_DA_forecasts_noIC, c("horizon", "model_id", "phen"), function(x) {
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
forecast_horizon_avg_noIC$model_id <- factor(forecast_horizon_avg_noIC$model_id, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))

#------------------------------------------------------------------------------------------------#
# Forecasts with IC on vs off
cb_friendly_2 <- c("#8C510A", "#BF812D", "#C7EAE5", "#35978F") #"#DFC27D", "#DEDEDE", "#8fd5cb"

#read in all forecasts with IC on
score_dir_yesIC <- arrow::SubTreeFileSystem$create(file.path(lake_directory,"scores/all_UC"))
all_DA_forecasts_yesIC <- arrow::open_dataset(score_dir_yesIC) |> 
  filter(!is.na(observation), variable == "temperature", horizon > 0.3) |>
  collect()    

#need to round horizon because they are in decimal form for all Jan 1 ref date
all_DA_forecasts_yesIC$horizon <- floor(all_DA_forecasts_yesIC$horizon)

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

#add new column for binning different layers of the water column
all_DA_forecasts_yesIC$layer <- ifelse(
  all_DA_forecasts_yesIC$depth <= 2.5, "surface",
  ifelse(all_DA_forecasts_yesIC$depth >= 8.5, "bottom", "middle"))

#forecast skill for each depth and horizon 
forecast_skill_layer_horizon_yesIC <-  plyr::ddply(all_DA_forecasts_yesIC, c("layer","horizon","phen", "model_id"), function(x) {
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
forecast_skill_layer_horizon_yesIC$model_id <- 
  factor(forecast_skill_layer_horizon_yesIC$model_id, 
         levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))

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

#add IC y/n column
forecast_skill_layer_horizon_noIC$IC <- "no"
forecast_skill_layer_horizon_yesIC$IC <- "yes"

forecast_horizon_avg_noIC$IC <- "n"
forecast_horizon_avg_yesIC$IC <- "y"

#fig s3
IC_off_rmse <- ggplot(subset(forecast_skill_layer_horizon_noIC, horizon > 0) ,
                     aes(horizon, RMSE, color=as.factor(model_id), aes=as.factor(IC))) +  
                ylab(expression("RMSE ("*degree*C*")")) + xlab("Horizon (days)")+
                geom_line() + theme_bw() + guides(color=guide_legend(title="DA frequency")) + 
                geom_hline(yintercept = 2, linetype="dotted") + ylim(0,3) +
                theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), 
                      legend.position = "right", legend.background = element_blank(),
                      panel.grid.minor = element_blank(), legend.key.size = unit(0.5, "lines"), 
                      panel.grid.major = element_blank(),legend.box.margin=margin(-10,-1,-10,-10),
                      legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"), 
                      legend.title = element_text(size=6), axis.text.x = element_text(vjust = 0.5,size=6), 
                      axis.text.y = element_text(size=6)) + #geom_point() +
                facet_grid(layer~phen, scales="free") + 
                scale_color_manual(values=cb_friendly_2)  

tag_facet2(IC_off_rmse, fontface = 1, hjust=0, size=3, tag_pool = c("a","b","c","d","e","f"))
#ggsave(file.path(lake_directory,"analysis/figures/UC_RMSEvsDAfreq_depth_facets_IC_off.jpg"),width=3.5, height=4)

#fig s4
IC_off_var <- ggplot(subset(forecast_skill_layer_horizon_noIC, horizon > 0) ,
                      aes(horizon, variance, color=as.factor(model_id), aes=as.factor(IC))) +  
                ylab("Variance") + xlab("Horizon (days)")+ ylim(0,8) +
                geom_line() + theme_bw() + guides(color=guide_legend(title="DA frequency")) + 
                theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), 
                      legend.position = "right", legend.background = element_blank(),
                      panel.grid.minor = element_blank(), legend.key.size = unit(0.5, "lines"), 
                      panel.grid.major = element_blank(),legend.box.margin=margin(-10,-1,-10,-10),
                      legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"), 
                      legend.title = element_text(size=6), axis.text.x = element_text(vjust = 0.5,size=6), 
                      axis.text.y = element_text(size=6)) + #geom_point() +
                facet_grid(layer~phen, scales="free") + 
                scale_color_manual(values=cb_friendly_2)  

tag_facet2(IC_off_var, fontface = 1, size=3, tag_pool = c("a","b","c","d","e","f"))  
#ggsave(file.path(lake_directory,"analysis/figures/UC_variancevsDAfreq_depth_facets_IC_off.jpg"),width=3.5, height=4)

#-------------------------------------------------------------------------------#
#FIGURE 9 - proportion of IC uncertainty vs horizon

#new dataframe to calculate the proportion of IC uncertainty across depths, horizons, model_id, and phen
UC_prop <- data.frame("layer" = forecast_skill_layer_horizon_yesIC$layer, 
                      "horizon" = forecast_skill_layer_horizon_yesIC$horizon, 
                      "phen" = forecast_skill_layer_horizon_yesIC$phen, 
                      "model_id" = forecast_skill_layer_horizon_yesIC$model_id,
                      "var" = forecast_skill_layer_horizon_yesIC$variance, 
                      "sd" = forecast_skill_layer_horizon_yesIC$sd)

UC_prop$prop_var <- ((forecast_skill_layer_horizon_yesIC$variance - 
                        forecast_skill_layer_horizon_noIC$variance) / 
                       forecast_skill_layer_horizon_yesIC$variance)

#set negative proportions to 0
UC_prop$prop_var[UC_prop$prop_var<0] <- 0

#order layers for plotting
UC_prop$layer <- factor(UC_prop$layer, levels=c("surface", "middle", "bottom"))
                           
fig9 <- ggplot(subset(UC_prop, horizon > 0) ,
               aes(horizon, prop_var, color=as.factor(model_id))) +  
  ylab("Proportion of total uncertainty contributed by initial conditions") + 
  geom_line() + theme_bw() +# geom_point() +
  guides(color=guide_legend(title="DA frequency")) + 
  xlab("Horizon (Days)")+ ylim(0, 0.85) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"),
        legend.position = "right",
        legend.background = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.5, "lines"), 
        panel.grid.major = element_blank(),
        legend.box.margin=margin(-10,-1,-10,-10),
        legend.text  = element_text(size = 6), 
        panel.spacing=unit(0, "cm"), 
        legend.title = element_text(size=6),
        axis.text.x = element_text(vjust = 0.5,size=6), 
        axis.text.y = element_text(size=6)) + 
  facet_grid(layer~phen, scales="free") + 
  scale_color_manual(values=cb_friendly_2) #+geom_point

tag_facet2(fig9, fontface = 1, size=3, tag_pool = c("a","b","c","d","e","f"))  
#ggsave(file.path(lake_directory,"analysis/figures/UC_prop_ICdvshorizon_depth_facets_var.jpg"),width=3.5, height=4)

mean(UC_prop$prop_var[UC_prop$horizon==1 & UC_prop$model_id=="Daily"])

mean(UC_prop$prop_var[UC_prop$horizon==1 & UC_prop$model_id=="Weekly"])
mean(UC_prop$prop_var[UC_prop$horizon==1 & UC_prop$model_id=="Fortnightly"])
mean(UC_prop$prop_var[UC_prop$horizon==1 & UC_prop$model_id=="Monthly"])

mean(UC_prop$prop_var[UC_prop$phen=="Mixed" & UC_prop$horizon>=10])
mean(UC_prop$prop_var[UC_prop$phen=="Stratified" & 
                        UC_prop$layer=="middle" &
                        UC_prop$horizon>=10])

mean(UC_prop$prop_var[UC_prop$phen=="Stratified" & UC_prop$layer=="middle" &
                        UC_prop$horizon>=10 & UC_prop$horizon<=20])

mean(UC_prop$prop_var[UC_prop$phen=="Mixed" & UC_prop$layer=="bottom" &
                        UC_prop$horizon==14 & UC_prop$model_id=="Monthly"])

#-----------------------------------------------------------------------------#
#tuned vs constant bad param runs - SI fig S5
  
  #read in all forecasts from tuned bad par run
  scores_tuned_bad_pars <- arrow::SubTreeFileSystem$create(file.path(lake_directory,"scores/tuned_bad_pars"))
  forecasts_tuned_bad_pars <- arrow::open_dataset(scores_tuned_bad_pars) |> 
    filter(!is.na(observation), variable == "temperature",horizon > 0.3) |>
    collect()    
  
  #need to round horizon because they are in decimal form for all Jan 1 ref date
  forecasts_tuned_bad_pars$horizon <- floor(forecasts_tuned_bad_pars$horizon)
  
  strat_date<- "2021-11-07"
  
  #add stratified vs mixed col
  forecasts_tuned_bad_pars$phen <- ifelse(forecasts_tuned_bad_pars$datetime <= as.POSIXct(strat_date) & 
                                            forecasts_tuned_bad_pars$datetime >="2021-03-12","Stratified", "Mixed")
  
  #remove n=6 days with ice-cover 
  forecasts_tuned_bad_pars <- forecasts_tuned_bad_pars[!(as.Date(forecasts_tuned_bad_pars$datetime) %in% c(as.Date("2021-01-10"), as.Date("2021-01-11"),as.Date("2021-01-30"),
                                                                                                  as.Date("2021-02-13"),as.Date("2021-02-14"),as.Date("2021-02-15"))),]
  #drop 11m completely because some rows were NA when water level was low
  forecasts_tuned_bad_pars <- forecasts_tuned_bad_pars[!(forecasts_tuned_bad_pars$depth==11),]
  
  #change model_id to be all uppercase
  forecasts_tuned_bad_pars$model_id <- str_to_title(forecasts_tuned_bad_pars$model_id)
  
  #only keep 2021 data
  forecasts_tuned_bad_pars <- forecasts_tuned_bad_pars[forecasts_tuned_bad_pars$datetime<="2021-12-31",]
  
  #add new column for binning different layers of the water column
  forecasts_tuned_bad_pars$layer <- ifelse(
    forecasts_tuned_bad_pars$depth <= 2.5, "surface",
    ifelse(forecasts_tuned_bad_pars$depth >= 8.5, "bottom", "middle"))
  
  #forecast skill for each depth and horizon
  forecasts_tuned_bad_pars_layer_horizon <-  plyr::ddply(forecasts_tuned_bad_pars, c("layer","horizon","phen", "model_id"), function(x) { #phen instead of datetime?
    data.frame(
      RMSE = sqrt(mean((x$mean - x$observation)^2, na.rm = TRUE)),
      MAE = mean(abs(x$mean - x$observation), na.rm = TRUE),
      pbias = 100 * (sum(x$mean - x$observation, na.rm = TRUE) / sum(x$observation, na.rm = TRUE)),
      CRPS = verification::crps(x$observation, as.matrix(x[, c(7,9)]))$CRPS,
      variance = (mean(x$sd))^2,
      sd = mean(x$sd)
    )
  }, .progress = plyr::progress_text(), .parallel = FALSE) 
  
  #rename daily to daily rerun
  forecasts_tuned_bad_pars_layer_horizon$model_id <- ifelse(
    forecasts_tuned_bad_pars_layer_horizon$model_id=="Daily_no_pars","tuned", 
    forecasts_tuned_bad_pars_layer_horizon$model_id)
  
  #read in all forecasts from constant bad parameter run
  scores_const_bad_pars <- arrow::SubTreeFileSystem$create(file.path(lake_directory,"scores/constant_bad_pars"))
  forecasts_const_bad_pars <- arrow::open_dataset(scores_const_bad_pars) |> 
    filter(!is.na(observation), variable == "temperature",horizon > 0.3) |>
    collect()    
  
  #need to round horizon because they are in decimal form for all Jan 1 ref date
  forecasts_const_bad_pars$horizon <- floor(forecasts_const_bad_pars$horizon)
  
  strat_date<- "2021-11-07"
  
  #add stratified vs mixed col
  forecasts_const_bad_pars$phen <- ifelse(forecasts_const_bad_pars$datetime <= as.POSIXct(strat_date) & 
                                            forecasts_const_bad_pars$datetime >="2021-03-12","Stratified", "Mixed")
  
  #remove n=6 days with ice-cover 
  forecasts_const_bad_pars <- forecasts_const_bad_pars[!(as.Date(forecasts_const_bad_pars$datetime) %in% c(as.Date("2021-01-10"), as.Date("2021-01-11"),as.Date("2021-01-30"),
                                                                                                           as.Date("2021-02-13"),as.Date("2021-02-14"),as.Date("2021-02-15"))),]
  #drop 11m completely because some rows were NA when water level was low
  forecasts_const_bad_pars <- forecasts_const_bad_pars[!(forecasts_const_bad_pars$depth==11),]
  
  #change model_id to be all uppercase
  forecasts_const_bad_pars$model_id <- str_to_title(forecasts_const_bad_pars$model_id)
  
  #only keep 2021 data
  forecasts_const_bad_pars <- forecasts_const_bad_pars[forecasts_const_bad_pars$datetime<="2021-12-31",]
  
  #add new column for binning different layers of the water column
  forecasts_const_bad_pars$layer <- ifelse(
    forecasts_const_bad_pars$depth <= 2.5, "surface",
    ifelse(forecasts_const_bad_pars$depth >= 8.5, "bottom", "middle"))
  
  #forecast skill for each depth and horizon
  forecasts_const_bad_pars_layer_horizon <-  plyr::ddply(forecasts_const_bad_pars, c("layer","horizon","phen", "model_id"), function(x) { #phen instead of datetime?
    data.frame(
      RMSE = sqrt(mean((x$mean - x$observation)^2, na.rm = TRUE)),
      MAE = mean(abs(x$mean - x$observation), na.rm = TRUE),
      pbias = 100 * (sum(x$mean - x$observation, na.rm = TRUE) / sum(x$observation, na.rm = TRUE)),
      CRPS = verification::crps(x$observation, as.matrix(x[, c(7,9)]))$CRPS,
      variance = (mean(x$sd))^2,
      sd = mean(x$sd)
    )
  }, .progress = plyr::progress_text(), .parallel = FALSE) 
  
  #rename daily to daily rerun
  forecasts_const_bad_pars_layer_horizon$model_id <- ifelse(
    forecasts_const_bad_pars_layer_horizon$model_id=="Daily_no_pars","constant", 
    forecasts_const_bad_pars_layer_horizon$model_id)
  
  #combine constant param and tuned param dfs
  daily_bad_pars <- dplyr::bind_rows(forecasts_tuned_bad_pars_layer_horizon,
                                              forecasts_const_bad_pars_layer_horizon)
  
  #order layer
  daily_bad_pars$layer <- factor(daily_bad_pars$layer, 
                                 levels=c( "surface", "middle", "bottom"))
  
  bad_par_daily <- ggplot(subset(daily_bad_pars, horizon > 0 &
                                model_id %in% c("constant", "tuned")) ,
                       aes(horizon, RMSE, color=as.factor(model_id))) +  
    ylab(expression("RMSE ("*degree*C*")")) + xlab("Horizon (days)")+
    geom_line() + theme_bw() + guides(color=guide_legend(title="parameters")) + 
    geom_hline(yintercept = 2, linetype="dotted") + #ylim(0,3.2) +
    theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), 
          legend.position = "right", legend.background = element_blank(),
          panel.grid.minor = element_blank(), legend.key.size = unit(0.5, "lines"), 
          panel.grid.major = element_blank(),legend.box.margin=margin(-10,-1,-10,-10),
          legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"), 
          legend.title = element_text(size=6), axis.text.x = element_text(vjust = 0.5,size=6), 
          axis.text.y = element_text(size=6)) + #geom_point() +
    facet_grid(layer~phen, scales="free") + 
    scale_color_manual(values=c("#CBCE91FF","#EA738DFF"))  
  
  tag_facet2(bad_par_daily, fontface = 1, size=3, tag_pool = c("a","b","c","d","e","f"))  
  #ggsave(file.path(lake_directory,"analysis/figures/Bad_pars_tuned_vs_constant_depth_facets.jpg"),width=3.5, height=4)
  
  mean(daily_bad_pars$RMSE[daily_bad_pars$model_id=="tuned"])
  mean(daily_bad_pars$sd[daily_bad_pars$model_id=="tuned"])
  
  mean(daily_bad_pars$RMSE[daily_bad_pars$model_id=="constant"])
  mean(daily_bad_pars$sd[daily_bad_pars$model_id=="constant"])
  
