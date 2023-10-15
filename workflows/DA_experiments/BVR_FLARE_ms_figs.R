#Figures for BVR FLARE ms
#09 Sep 2022 HLW
library(tidyverse)
library(lubridate)
library(ggnewscale)
library(viridis)
library(padr)

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
forecast_site <- "bvre"
configure_run_file <- "configure_run.yml"
config_files <- "configure_flare.yml"
config_set_name <- "DA_experiments"
setwd(lake_directory)

#read in all forecasts 
score_dir <- arrow::SubTreeFileSystem$create(file.path(lake_directory,"scores/all_UC"))
all_DA_forecasts <- arrow::open_dataset(score_dir) |> 
  filter(!is.na(observation), variable == "temperature", horizon >= 0) |>
         collect() |>   
  filter(as.Date(reference_datetime) > "2020-12-31") 

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

#add new column for binning different layers of the water column
all_DA_forecasts$layer <- ifelse(
  all_DA_forecasts$depth <= 2.5, "surface",
  ifelse(all_DA_forecasts$depth >= 8.5, "bottom", "middle"))

#now summarise by layer and horizon
forecast_skill_layer_horizon <-  plyr::ddply(
  all_DA_forecasts, c("layer","horizon", "phen", "model_id"), 
  function(x) {
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
forecast_skill_layer_horizon$model_id <- 
  factor(forecast_skill_layer_horizon$model_id, levels=c(
    "Daily", "Weekly", "Fortnightly", "Monthly"))

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

#combine constant param and tuned param dfs
#test <- dplyr::bind_rows(forecasts_tuned_bad_pars_depth_horizon_avg,
#                         forecast_avg)

#Figure to answer Q1 - aggregated RMSE across depths and seasons
ggplot(subset(forecast_avg, horizon > 0 & #test
                model_id %in% c("Daily","Weekly","Fortnightly","Monthly")) , #tuned
       aes(horizon, RMSE, color=as.factor(model_id))) +  
  ylab(expression("RMSE ("*degree*C*")")) + xlab("Horizon (days)")+
  geom_line() + theme_bw() + guides(color=guide_legend(title="DA frequency")) + 
  geom_hline(yintercept = 2, linetype="dotted") + ylim(0,2.2) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = "right",
        legend.background = element_blank(),legend.direction = "vertical", panel.grid.minor = element_blank(),
        legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),legend.box.margin=margin(-10,-1,-10,-10),
        legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"), legend.title = element_text(size=6),
        axis.text.x = element_text(vjust = 0.5,size=6), axis.text.y = element_text(size=6)) +
  scale_color_manual(values=cb_friendly_2)# +geom_point()
#ggsave(file.path(lake_directory,"analysis/figures/RMSEvshorizon_fig6_Q1_daily_badpars.jpg"),width=4, height=3)

#reorder layers
forecast_skill_layer_horizon$layer <- 
  factor(forecast_skill_layer_horizon$layer, levels = c("surface", "middle","bottom"))

#FIGURE 7 - RMSE vs horizon
fig7 <- ggplot(subset(forecast_skill_layer_horizon, horizon > 0 & #test
                        model_id %in% c("Daily","Weekly","Fortnightly","Monthly")) , #tuned
               aes(horizon, RMSE, color=as.factor(model_id))) +  
  ylab(expression("RMSE ("*degree*C*")")) + xlab("Horizon (days)")+
   geom_line() + theme_bw() + guides(color=guide_legend(title="DA frequency")) + 
  geom_hline(yintercept = 2, linetype="dotted") + ylim(0,3.2) +
  theme(text = element_text(size=8), 
        axis.text = element_text(size=6, color="black"), 
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
        axis.text.y = element_text(size=6)) + #geom_point() +
  facet_grid(layer~phen, scales="free") + 
  scale_color_manual(values=cb_friendly_2) 

tag_facet2(fig7, fontface = 1, size=3,
           tag_pool = c("a","b","c","d","e","f"))
#ggsave(file.path(lake_directory,"analysis/figures/RMSEvshorizon_depth_facets_fig7.jpg"),width=3.5, height=4)


fig8 <- ggplot(subset(forecast_skill_layer_horizon, horizon > 0) ,
               aes(horizon, variance, color=as.factor(model_id))) +  
  ylab("Total variance") + xlab("Horizon (days)")+
  geom_line() + theme_bw() + guides(color=guide_legend(title="DA frequency")) + 
  ylim(0,8) + 
  theme(text = element_text(size=8), 
        axis.text = element_text(size=6, color="black"), 
        legend.position = "right",
        legend.background = element_blank(),
        legend.direction = "vertical", 
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
  scale_color_manual(values=cb_friendly_2) 

tag_facet2(fig8, fontface = 1, size=3,
           tag_pool = c("a","b","c","d","e","f"))
#ggsave(file.path(lake_directory,"analysis/figures/Varvshorizon_depth_facets_fig8.jpg"),width=3.5, height=4)

range(forecast_skill_layer_horizon$variance[
  forecast_skill_layer_horizon$horizon<=7 &
        forecast_skill_layer_horizon$layer=="surface"])

range(forecast_skill_layer_horizon$variance[
  forecast_skill_layer_horizon$horizon<=7 &
    forecast_skill_layer_horizon$layer=="middle"])

range(forecast_skill_layer_horizon$variance[
  forecast_skill_layer_horizon$horizon<=7 &
    forecast_skill_layer_horizon$layer=="bottom"])

mean(forecast_skill_layer_horizon$RMSE)
mean(forecast_skill_layer_horizon$sd)

mean(forecast_skill_layer_horizon$variance[forecast_skill_layer_horizon$horizon==1 &
                                             forecast_skill_layer_horizon$model_id=="Daily"])

mean(forecast_skill_layer_horizon$variance[forecast_skill_layer_horizon$horizon==1 &
                                             forecast_skill_layer_horizon$model_id=="Monthly"])

mean(forecast_skill_layer_horizon$RMSE[forecast_skill_layer_horizon$horizon>=8 & 
                                         forecast_skill_layer_horizon$model_id=="Daily" &
                                         forecast_skill_layer_horizon$layer=="bottom"])

mean(forecast_skill_layer_horizon$RMSE[forecast_skill_layer_horizon$horizon>=8 & 
                                         forecast_skill_layer_horizon$model_id=="Weekly"&
                                         forecast_skill_layer_horizon$layer=="bottom"])

mean(forecast_skill_layer_horizon$RMSE[forecast_skill_layer_horizon$horizon==1])
mean(forecast_skill_layer_horizon$sd[forecast_skill_layer_horizon$horizon==1])

mean(forecast_skill_layer_horizon$RMSE[forecast_skill_layer_horizon$horizon==7])
mean(forecast_skill_layer_horizon$sd[forecast_skill_layer_horizon$horizon==7])

mean(forecast_skill_layer_horizon$RMSE[forecast_skill_layer_horizon$horizon==35])
mean(forecast_skill_layer_horizon$sd[forecast_skill_layer_horizon$horizon==35])

mean(forecast_skill_layer_horizon$RMSE[forecast_skill_layer_horizon$phen=="Mixed"])
mean(forecast_skill_layer_horizon$sd[forecast_skill_layer_horizon$phen=="Mixed"])

mean(forecast_skill_layer_horizon$RMSE[forecast_skill_layer_horizon$phen=="Stratified"])
mean(forecast_skill_layer_horizon$sd[forecast_skill_layer_horizon$phen=="Stratified"])

mean(forecast_skill_layer_horizon$RMSE[forecast_skill_layer_horizon$layer=="surface"])
mean(forecast_skill_layer_horizon$sd[forecast_skill_layer_horizon$layer=="surface"])

mean(forecast_skill_layer_horizon$RMSE[forecast_skill_layer_horizon$layer=="middle"])
mean(forecast_skill_layer_horizon$sd[forecast_skill_layer_horizon$layer=="middle"])

mean(forecast_skill_layer_horizon$RMSE[forecast_skill_layer_horizon$layer=="bottom"])
mean(forecast_skill_layer_horizon$sd[forecast_skill_layer_horizon$layer=="bottom"])

max(forecast_skill_layer_horizon$RMSE[forecast_skill_layer_horizon$layer=="bottom" & forecast_skill_layer_horizon$phen=="Stratified"])

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

#change first 10 to 9m for 6jul-13Jul bc of rounding issues with changing water level - super hacky...
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

layers <- c('surface','surface','surface',
            'middle','middle','middle','middle',
            'bottom','bottom','bottom','bottom')
names(layers) <- unique(all_da_sub_final$depth)

fig4_poc <- ggplot(subset(all_da_sub_final, depth %in% c(1,5,9)), aes(horizon, forecast_mean, color=model_id)) + 
  geom_ribbon(data=subset(all_da_sub_final, depth %in% c(1,5,9) & model_id=="Monthly"), 
              aes(x=date, y = forecast_mean, ymin = forecast_mean-forecast_sd, 
                  ymax = forecast_mean+forecast_sd), color=cb_friendly_2[4], 
              fill=cb_friendly_2[4], alpha=0.4) +
  geom_ribbon(data=subset(all_da_sub_final, depth %in% c(1,5,9) & model_id=="Fortnightly"), 
              aes(x=date, y = forecast_mean, ymin = forecast_mean-forecast_sd, 
                  ymax = forecast_mean+forecast_sd), color=cb_friendly_2[3], 
              fill=cb_friendly_2[3], alpha=0.4) +
  geom_ribbon(data=subset(all_da_sub_final, depth %in% c(1,5,9) & model_id=="Weekly"), 
              aes(x=date, y = forecast_mean, ymin = forecast_mean-forecast_sd, 
                  ymax = forecast_mean+forecast_sd), color=cb_friendly_2[2], 
              fill=cb_friendly_2[2], alpha=0.4)  +
  geom_ribbon(data=subset(all_da_sub_final, depth %in% c(1,5,9) & model_id=="Daily"), 
              aes(x=date, y = forecast_mean, ymin = forecast_mean-forecast_sd, 
                  ymax = forecast_mean+forecast_sd), color=cb_friendly_2[1], 
              fill=cb_friendly_2[1], alpha=0.4) + theme_bw() + 
  scale_fill_manual(values=rev(cb_friendly_2), 
                    labels=c("Daily","Weekly","Fortnightly","Monthly")) + 
  guides(fill="none") + new_scale_fill() +
  geom_point(data=subset(obs, depth %in% c(1,5,9)), 
             aes(x=datetime, y=observation, size=siz, fill=col), 
             pch=21, color="black") + 
  scale_fill_manual(values=c(rev(cb_friendly_2)), 
                    labels=c("Daily","Weekly","Fortnightly","Monthly", "")) +
  scale_x_date(date_breaks = "1 week", date_labels = "%b %d", 
               limits = c(as.Date("2021-06-25"),NA)) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), 
        legend.position = "right", legend.box.margin=margin(-10,-1,-10,-10),
        legend.background = element_blank(), panel.grid.minor = element_blank(), 
        legend.key=element_rect(fill=NA), plot.margin = unit(c(0,0.05,-0.2,0), "cm"),
        legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 4.5),legend.text  = element_text(size = 4.5), 
        panel.spacing=unit(0, "cm"), axis.text.y = element_text(size=6),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6)) +
  facet_wrap(~depth, ncol = 1, scales = "free_y",labeller = labeller(depth = layers)) +
  ylab(expression("Water temperature ("*degree*C*")")) + xlab("")  + 
  scale_color_manual(values=rev(cb_friendly_2)) + 
  scale_size_manual(values=c(0.8,0.1)) +
  #geom_text(x=as.Date("2021-07-18"), y=24, label="Past", color="black", size=2) + 
  #geom_text(x=as.Date("2021-07-26"), y=24, label="Future", color="black", size=2) +
  guides(fill = guide_legend(title="DA frequency", 
                             override.aes = list(alpha=1,
                                                pch=c(rep(21,4),NA),
                                                fill=c(cb_friendly_2,NA)),reverse = FALSE), 
         color="none", size="none")
tag_facet2(fig4_poc, fontface = 1, hjust=0.7, size=3,
           tag_pool = c("a","b","c"),  x = as.Date("2021-06-24"))  
#ggsave(file.path(lake_directory,"analysis/figures/forecast_ProofOfConcept_fig4.jpg"), width=3.5, height=4) 


#% uncertainty at 0-day horizon for daily and monthly forecasts
mean(all_da_sub_final$forecast_upper_95[all_da_sub_final$horizon==1 & all_da_sub_final$model_id=="Daily"]) -
mean(all_da_sub_final$forecast_lower_95[all_da_sub_final$horizon==1 & all_da_sub_final$model_id=="Daily"])

mean(all_da_sub_final$forecast_upper_95[all_da_sub_final$horizon==29 & all_da_sub_final$model_id=="Monthly"]) -
  mean(all_da_sub_final$forecast_lower_95[all_da_sub_final$horizon==29 & all_da_sub_final$model_id=="Monthly"])


#uncertainty range across depths and mixed vs stratified periods
((mean(all_DA_forecasts$quantile97.5[all_DA_forecasts$phen=="Mixed" & all_DA_forecasts$model_id=="Monthly"]) -
    mean(all_DA_forecasts$quantile02.5[all_DA_forecasts$phen=="Mixed" & all_DA_forecasts$model_id=="Monthly"])) -
    
    (mean(all_DA_forecasts$quantile97.5[all_DA_forecasts$phen=="Mixed" & all_DA_forecasts$model_id=="Daily"]) -
       mean(all_DA_forecasts$quantile02.5[all_DA_forecasts$phen=="Mixed" & all_DA_forecasts$model_id=="Daily"]))) /
  
  mean((mean(all_DA_forecasts$quantile97.5[all_DA_forecasts$phen=="Mixed" & all_DA_forecasts$model_id=="Monthly"]) -
          mean(all_DA_forecasts$quantile02.5[all_DA_forecasts$phen=="Mixed" & all_DA_forecasts$model_id=="Monthly"])),
       (mean(all_DA_forecasts$quantile97.5[all_DA_forecasts$phen=="Mixed" & all_DA_forecasts$model_id=="Daily"]) -
          mean(all_DA_forecasts$quantile02.5[all_DA_forecasts$phen=="Mixed" & all_DA_forecasts$model_id=="Daily"]))) *100


((mean(all_DA_forecasts$quantile97.5[all_DA_forecasts$phen=="Stratified" & all_DA_forecasts$model_id=="Monthly"]) -
    mean(all_DA_forecasts$quantile02.5[all_DA_forecasts$phen=="Stratified" & all_DA_forecasts$model_id=="Monthly"])) -
    
    (mean(all_DA_forecasts$quantile97.5[all_DA_forecasts$phen=="Stratified" & all_DA_forecasts$model_id=="Daily"]) -
       mean(all_DA_forecasts$quantile02.5[all_DA_forecasts$phen=="Stratified" & all_DA_forecasts$model_id=="Daily"]))) / 
  
  mean((mean(all_DA_forecasts$quantile97.5[all_DA_forecasts$phen=="Stratified" & all_DA_forecasts$model_id=="Monthly"]) -
          mean(all_DA_forecasts$quantile02.5[all_DA_forecasts$phen=="Stratified" & all_DA_forecasts$model_id=="Monthly"])),
       (mean(all_DA_forecasts$quantile97.5[all_DA_forecasts$phen=="Stratified" & all_DA_forecasts$model_id=="Daily"]) -
          mean(all_DA_forecasts$quantile02.5[all_DA_forecasts$phen=="Stratified" & all_DA_forecasts$model_id=="Daily"]))) *100

#--------------------------------------------------------------------------------#
#Fig 3 - water temp phenology fig
# Read in FLARE observations ----
config <- FLAREr::set_configuration(configure_run_file, lake_directory, 
                                    config_set_name = config_set_name, 
                                    sim_name = NA)

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

ggplot(sub) +   theme_bw() +
  geom_rect(data = phen, aes(fill = "Mixed"), xmin=-Inf ,xmax = as.Date("2021-03-11"),
            ymin = -Inf, ymax = Inf, inherit.aes = FALSE) + 
  geom_rect(data = phen, aes(fill = "Stratified"), xmin=as.Date("2021-03-12") ,
            xmax = as.Date("2021-11-07"), ymin = -Inf, ymax = Inf, inherit.aes = FALSE)+
  geom_rect(data = phen, aes(fill = "Mixed"), xmin=as.Date("2021-11-08") ,
            xmax = Inf, ymin = -Inf, ymax = Inf, inherit.aes = FALSE) +
  geom_line(aes(Date, as.numeric(observation), color = factor(depth)), 
            linewidth=0.4) + ylab(expression("Water temperature ("*degree*C*")")) + 
  xlab("") +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), 
        legend.position = "right", legend.background = element_blank(),
        legend.key = element_blank(), legend.key.height = unit(0.3,"cm"), 
        legend.key.width = unit(0.4,"cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_fill_manual('', values = c('gray','white')) +
  scale_color_viridis(option="B", discrete="TRUE", direction=-1) +
  guides(color=guide_legend("Depth (m)"),
         fill= guide_legend("Period",order = 1, 
                            override.aes= list(color="black", lwd=0.1)))
#ggsave(file.path(lake_directory,"analysis/figures/2021_watertemp_mixedVstratified_fig3.jpg"), width=4, height=3)

range(sub$observation[sub$datetime>= "2021-03-12" & sub$datetime<= "2021-11-07"])

#-------------------------------------------------------------------------------#
#parameter evolution figs 
#read in all forecasts 
params_dir <- arrow::SubTreeFileSystem$create(file.path(lake_directory,"scores/all_UC"))
params <- arrow::open_dataset(params_dir) |>  
  filter(variable %in% c("lw_factor","zone1temp","zone2temp"), horizon >=0) |>  
  collect() 

#change datetime format
params$datetime <- as.Date(params$datetime)

#change model_id to be all uppercase
params$model_id <- str_to_title(params$model_id)

#change DA factor order
params$model_id <- factor(params$model_id, levels = c("Daily", "Weekly","Fortnightly","Monthly"))

variable <- c("longwave","hypo_sed_temp","epi_sed_temp")
names(variable) <- c("lw_factor","zone1temp","zone2temp")

fig5 <- ggplot(subset(params,horizon==1), aes(as.Date(datetime), mean, color=model_id)) + theme_bw() +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position ="right",
        legend.background = element_blank(), panel.grid.minor = element_blank(), legend.box.margin=margin(-10,-1,-10,-10),
        plot.margin = unit(c(0,0.05,-0.2,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 6),legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  scale_color_manual(values=cb_friendly_2) +
  scale_fill_manual(values=cb_friendly_2) +
  facet_wrap(~variable, scales="free_y",ncol=1, labeller = labeller(variable = variable)) +
  ylab("Parameter value")+ xlab("")+
  scale_x_date(date_labels = "%b") + #ylab(expression("Temperature ("*~degree*C*")")) 
  geom_ribbon(aes(y = mean, ymin = mean-sd, ymax = mean+sd, color=model_id, fill=model_id), alpha=0.5) +
  guides(fill = guide_legend(title="DA frequency", override.aes = list(alpha=1)), color="none")

tag_facet2(fig5, fontface = 1, x=as.Date("2021-01-01"), hjust=0.7, size=3,
           tag_pool = c("a","b","c"))
#ggsave(file.path(lake_directory,"analysis/figures/paramevolvstime_1day_fig5.jpg"),width=3.5, height=4)

#figuring out the date that DA parameters diverge
mean(params$mean[params$variable=="lw_factor" & params$model_id=="Daily" & params$datetime >= "2021-04-01"])

mean(params$mean[params$variable=="lw_factor" & params$model_id=="Weekly" & params$datetime >= "2021-04-01"])
mean(params$mean[params$variable=="lw_factor" & params$model_id=="Fortnightly" & params$datetime >= "2021-04-01"])
mean(params$mean[params$variable=="lw_factor" & params$model_id=="Monthly" & params$datetime >= "2021-04-01"])

#ranges for daily vs non-daily sed temp params
range(params$mean[params$variable=="zone2temp" & params$model_id!="Daily" & params$datetime >= "2021-04-01"])
range(params$mean[params$variable=="zone2temp" & params$model_id=="Daily" & params$datetime >= "2021-04-01"])

range(params$mean[params$variable=="zone1temp" & params$model_id!="Daily" & params$datetime >= "2021-04-01"])
range(params$mean[params$variable=="zone1temp" & params$model_id=="Daily" & params$datetime >= "2021-04-01"])


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

#-------------------------------------------------------------------------------#
#SI figs

#CRPS fig across mixed vs strat, depths, and horizons (fig 4 but for CRPS)

#reorder layers
forecast_skill_layer_horizon$layer <- 
  factor(forecast_skill_layer_horizon$layer, levels = c("surface", "middle","bottom"))

figs2 <- ggplot(subset(forecast_skill_layer_horizon, horizon > 0) ,
               aes(horizon, CRPS, color=as.factor(model_id))) +  
  ylab(expression("CRPS ("*degree*C*")")) + xlab("Horizon (days)")+
  geom_line() + theme_bw() + guides(color=guide_legend(title="DA frequency")) + 
  ylim(0,2) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), 
        legend.position = "right", legend.background = element_blank(),
        panel.grid.minor = element_blank(), legend.key.size = unit(0.5, "lines"), 
        panel.grid.major = element_blank(),legend.box.margin=margin(-10,-1,-10,-10),
        legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"), 
        legend.title = element_text(size=6), axis.text.x = element_text(vjust = 0.5,size=6), 
        axis.text.y = element_text(size=6)) + #geom_point() +
  facet_grid(layer~phen, scales="free") + 
  scale_color_manual(values=cb_friendly_2) 

tag_facet2(figs2, fontface = 1, size=3,
           tag_pool = c("a","b","c","d","e","f"))

#ggsave(file.path(lake_directory,"analysis/figures/CRPSvsDAfreq_depth_facets_figs2.jpg"),width=3.5, height=4)

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
        axis.text.x = element_text(vjust = 0.5,size=6), 
        axis.text.y = element_text(size=6)) 
#ggsave(file.path(lake_directory,"analysis/figures/DA_freq_2021_S1.jpg"), width=4, height=3.5) 
