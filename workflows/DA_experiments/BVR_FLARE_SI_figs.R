#BVR FLARE SI figs

#load libraries
pacman::p_load(reshape2)

#inflation parameter figs (fixed, 1.02, 1.04)
#and different start days (27nov, 24nov, 22nov)

#read in all forecast csvs 
daily_27nov <- list.files(file.path(lake_directory,"analysis/summary_files/27nov_start/daily"), pattern="csv", full.names=TRUE)
daily_27nov <- lapply(daily_27nov, read_csv) %>% bind_rows() %>% mutate(DA = "Daily")

daily_24nov <- list.files(file.path(lake_directory,"analysis/summary_files/24nov_start/daily"), pattern="csv", full.names=TRUE)
daily_24nov <- lapply(daily_24nov, read_csv) %>% bind_rows() %>% mutate(DA = "Daily")

daily_22nov <- list.files(file.path(lake_directory,"analysis/summary_files/22nov_start/daily"), pattern="csv", full.names=TRUE)
daily_22nov <- lapply(daily_22nov, read_csv) %>% bind_rows() %>% mutate(DA = "Daily")

daily_fixed <- list.files(file.path(lake_directory,"analysis/summary_files/fixed_params/daily"), pattern="csv", full.names=TRUE)
daily_fixed <- lapply(daily_fixed, read_csv) %>% bind_rows() %>% mutate(DA = "Daily")

daily_1.04 <- list.files(file.path(lake_directory,"analysis/summary_files/inflat_1.04/daily"), pattern="csv", full.names=TRUE)
daily_1.04 <- lapply(daily_1.04, read_csv) %>% bind_rows() %>% mutate(DA = "Daily")

weekly_27nov <- list.files(file.path(lake_directory,"analysis/summary_files/27nov_start/weekly"), pattern="csv", full.names=TRUE)
weekly_27nov <- lapply(weekly_27nov, read_csv) %>% bind_rows() %>% mutate(DA = "Weekly")

weekly_24nov <- list.files(file.path(lake_directory,"analysis/summary_files/24nov_start/weekly"), pattern="csv", full.names=TRUE)
weekly_24nov <- lapply(weekly_24nov, read_csv) %>% bind_rows() %>% mutate(DA = "Weekly")

weekly_22nov <- list.files(file.path(lake_directory,"analysis/summary_files/22nov_start/weekly"), pattern="csv", full.names=TRUE)
weekly_22nov <- lapply(weekly_22nov, read_csv) %>% bind_rows() %>% mutate(DA = "Weekly")

weekly_fixed <- list.files(file.path(lake_directory,"analysis/summary_files/fixed_params/weekly"), pattern="csv", full.names=TRUE)
weekly_fixed <- lapply(weekly_fixed, read_csv) %>% bind_rows() %>% mutate(DA = "Weekly")

weekly_1.04 <- list.files(file.path(lake_directory,"analysis/summary_files/inflat_1.04/weekly"), pattern="csv", full.names=TRUE)
weekly_1.04 <- lapply(weekly_1.04, read_csv) %>% bind_rows() %>% mutate(DA = "Weekly")

fortnightly_27nov <- list.files(file.path(lake_directory,"analysis/summary_files/27nov_start/fortnightly"), pattern="csv", full.names=TRUE)
fortnightly_27nov <- lapply(fortnightly_27nov, read_csv) %>% bind_rows() %>% mutate(DA = "Fortnightly")

fortnightly_24nov <- list.files(file.path(lake_directory,"analysis/summary_files/24nov_start/fortnightly"), pattern="csv", full.names=TRUE)
fortnightly_24nov <- lapply(fortnightly_24nov, read_csv) %>% bind_rows() %>% mutate(DA = "Fortnightly")

fortnightly_22nov <- list.files(file.path(lake_directory,"analysis/summary_files/22nov_start/fortnightly"), pattern="csv", full.names=TRUE)
fortnightly_22nov <- lapply(fortnightly_22nov, read_csv) %>% bind_rows() %>% mutate(DA = "Fortnightly")

fortnightly_fixed <- list.files(file.path(lake_directory,"analysis/summary_files/fixed_params/fortnightly"), pattern="csv", full.names=TRUE)
fortnightly_fixed <- lapply(fortnightly_fixed, read_csv) %>% bind_rows() %>% mutate(DA = "Fortnightly")

fortnightly_1.04 <- list.files(file.path(lake_directory,"analysis/summary_files/inflat_1.04/fortnightly"), pattern="csv", full.names=TRUE)
fortnightly_1.04 <- lapply(fortnightly_1.04, read_csv) %>% bind_rows() %>% mutate(DA = "Fortnightly")

monthly_27nov <- list.files(file.path(lake_directory,"analysis/summary_files/27nov_start/monthly"), pattern="csv", full.names=TRUE)
monthly_27nov <- lapply(monthly_27nov, read_csv) %>% bind_rows() %>% mutate(DA = "Monthly")

monthly_24nov <- list.files(file.path(lake_directory,"analysis/summary_files/24nov_start/monthly"), pattern="csv", full.names=TRUE)
monthly_24nov <- lapply(monthly_24nov, read_csv) %>% bind_rows() %>% mutate(DA = "Monthly")

monthly_22nov <- list.files(file.path(lake_directory,"analysis/summary_files/22nov_start/monthly"), pattern="csv", full.names=TRUE)
monthly_22nov <- lapply(monthly_22nov, read_csv) %>% bind_rows() %>% mutate(DA = "Monthly")

monthly_fixed <- list.files(file.path(lake_directory,"analysis/summary_files/fixed_params/monthly"), pattern="csv", full.names=TRUE)
monthly_fixed <- lapply(monthly_fixed, read_csv) %>% bind_rows() %>% mutate(DA = "Monthly")

monthly_1.04 <- list.files(file.path(lake_directory,"analysis/summary_files/inflat_1.04/monthly"), pattern="csv", full.names=TRUE)
monthly_1.04 <- lapply(monthly_1.04, read_csv) %>% bind_rows() %>% mutate(DA = "Monthly")


detach(dplyr)
library(plyr)

#merge all DA frequencies
all_27nov <- rbind(daily_27nov, weekly_27nov, fortnightly_27nov, monthly_27nov)
all_24nov <- rbind(daily_24nov, weekly_24nov, fortnightly_24nov, monthly_24nov)
all_22nov <- rbind(daily_22nov, weekly_22nov, fortnightly_22nov, monthly_22nov)
all_fixed <- rbind(daily_fixed, weekly_fixed, fortnightly_fixed, monthly_fixed)
all_1.04 <- rbind(daily_1.04, weekly_1.04, fortnightly_1.04, monthly_1.04)


#round depths to nearest m
all_27nov$depth <- ceiling(all_27nov$depth)
all_24nov$depth <- ceiling(all_24nov$depth)
all_22nov$depth <- ceiling(all_22nov$depth)
all_fixed$depth <- ceiling(all_fixed$depth)
all_1.04$depth <- ceiling(all_1.04$depth)


#add date column --> forecast_date + horizon
all_27nov$date <- all_27nov$forecast_date + all_27nov$horizon
all_24nov$date <- all_24nov$forecast_date + all_24nov$horizon
all_22nov$date <- all_22nov$forecast_date + all_22nov$horizon
all_fixed$date <- all_fixed$forecast_date + all_fixed$horizon
all_1.04$date <-  all_1.04$forecast_date + all_1.04$horizon

strat_date<- "2021-11-07"

#add stratified vs mixed col
all_27nov$phen <- ifelse(all_27nov$date <= as.Date(strat_date) & 
                           all_27nov$date >="2021-03-13","Stratified", "Mixed")

all_24nov$phen <- ifelse(all_24nov$date <= as.Date(strat_date) & 
                           all_24nov$date >="2021-03-13","Stratified", "Mixed")

all_22nov$phen <- ifelse(all_22nov$date <= as.Date(strat_date) & 
                           all_22nov$date >="2021-03-13","Stratified", "Mixed")

all_fixed$phen <- ifelse(all_fixed$date <= as.Date(strat_date) & 
                           all_fixed$date >="2021-03-13","Stratified", "Mixed")

all_1.04$phen <- ifelse(all_1.04$date <= as.Date(strat_date) & 
                          all_1.04$date >="2021-03-13","Stratified", "Mixed")

#only keep 2021 data
all_27nov <- all_27nov[all_27nov$date<="2021-12-31",]
all_24nov <- all_24nov[all_24nov$date<="2021-12-31",]
all_22nov <- all_22nov[all_22nov$date<="2021-12-31",]
all_fixed <- all_fixed[all_fixed$date<="2021-12-31",]
all_1.04 <- all_1.04[all_1.04$date<="2021-12-31",]

#drop 11m completely because some rows were NA when water level was low
all_27nov <- all_27nov[!(all_27nov$depth==11),]
all_24nov <- all_24nov[!(all_24nov$depth==11),]
all_22nov <- all_22nov[!(all_22nov$depth==11),]
all_fixed <- all_fixed[!(all_fixed$depth==11),]
all_1.04 <- all_1.04[!(all_1.04$depth==11),]

#forecast skill for each depth and horizon
forecast_skill_depth_horizon_27nov <-  plyr::ddply(all_27nov, c("depth","horizon","phen", "DA"), function(x) {
  data.frame(
    RMSE = sqrt(mean((x$mean - x$obs)^2, na.rm = TRUE)),
    MAE = mean(abs(x$mean - x$obs), na.rm = TRUE),
    pbias = 100 * (sum(x$mean - x$obs, na.rm = TRUE) / sum(x$obs, na.rm = TRUE)),
    CRPS = verification::crps(x$obs, as.matrix(x[, 4:5]))$CRPS,
    variance = (mean(x$sd))^2
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

forecast_skill_depth_horizon_24nov <-  plyr::ddply(all_24nov, c("depth","horizon","phen", "DA"), function(x) {
  data.frame(
    RMSE = sqrt(mean((x$mean - x$obs)^2, na.rm = TRUE)),
    MAE = mean(abs(x$mean - x$obs), na.rm = TRUE),
    pbias = 100 * (sum(x$mean - x$obs, na.rm = TRUE) / sum(x$obs, na.rm = TRUE)),
    CRPS = verification::crps(x$obs, as.matrix(x[, 4:5]))$CRPS,
    variance = (mean(x$sd))^2
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

forecast_skill_depth_horizon_22nov <-  plyr::ddply(all_22nov, c("depth","horizon","phen", "DA"), function(x) {
  data.frame(
    RMSE = sqrt(mean((x$mean - x$obs)^2, na.rm = TRUE)),
    MAE = mean(abs(x$mean - x$obs), na.rm = TRUE),
    pbias = 100 * (sum(x$mean - x$obs, na.rm = TRUE) / sum(x$obs, na.rm = TRUE)),
    CRPS = verification::crps(x$obs, as.matrix(x[, 4:5]))$CRPS,
    variance = (mean(x$sd))^2
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

forecast_skill_depth_horizon_fixed <-  plyr::ddply(all_fixed, c("depth","horizon","phen", "DA"), function(x) {
  data.frame(
    RMSE = sqrt(mean((x$mean - x$obs)^2, na.rm = TRUE)),
    MAE = mean(abs(x$mean - x$obs), na.rm = TRUE),
    pbias = 100 * (sum(x$mean - x$obs, na.rm = TRUE) / sum(x$obs, na.rm = TRUE)),
    CRPS = verification::crps(x$obs, as.matrix(x[, 4:5]))$CRPS,
    variance = (mean(x$sd))^2
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

forecast_skill_depth_horizon_1.04 <-  plyr::ddply(all_1.04, c("depth","horizon","phen", "DA"), function(x) {
  data.frame(
    RMSE = sqrt(mean((x$mean - x$obs)^2, na.rm = TRUE)),
    MAE = mean(abs(x$mean - x$obs), na.rm = TRUE),
    pbias = 100 * (sum(x$mean - x$obs, na.rm = TRUE) / sum(x$obs, na.rm = TRUE)),
    CRPS = verification::crps(x$obs, as.matrix(x[, 4:5]))$CRPS,
    variance = (mean(x$sd))^2
  )
}, .progress = plyr::progress_text(), .parallel = FALSE) 

#order DA frequencies
forecast_skill_depth_horizon_27nov$DA <- factor(forecast_skill_depth_horizon_27nov$DA, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))
forecast_skill_depth_horizon_24nov$DA <- factor(forecast_skill_depth_horizon_24nov$DA, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))
forecast_skill_depth_horizon_22nov$DA <- factor(forecast_skill_depth_horizon_22nov$DA, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))
forecast_skill_depth_horizon_fixed$DA <- factor(forecast_skill_depth_horizon_fixed$DA, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))
forecast_skill_depth_horizon_1.04$DA <- factor(forecast_skill_depth_horizon_1.04$DA, levels=c("Daily", "Weekly", "Fortnightly", "Monthly"))


#------------------------------------------------------------------------------------------------#
# Start date forecast skill comparison

#add new column for start date / parameter
forecast_skill_depth_horizon_27nov$da_start_date <- "27nov"
forecast_skill_depth_horizon_24nov$da_start_date <- "24nov"
forecast_skill_depth_horizon_22nov$da_start_date <- "22nov"
forecast_skill_depth_horizon_fixed$da_start_date <- "27nov"
forecast_skill_depth_horizon_1.04$da_start_date <- "27nov"


forecast_skill_depth_horizon_27nov$param <- "1.02"
forecast_skill_depth_horizon_24nov$param <- "1.02"
forecast_skill_depth_horizon_22nov$param <- "1.02"
forecast_skill_depth_horizon_fixed$param <- "fixed"
forecast_skill_depth_horizon_1.04$param <- "1.04"

#combine to make a massive df
start_dates <- rbind(forecast_skill_depth_horizon_27nov,forecast_skill_depth_horizon_24nov,forecast_skill_depth_horizon_22nov)
params <- rbind(forecast_skill_depth_horizon_27nov,forecast_skill_depth_horizon_fixed,forecast_skill_depth_horizon_1.04)

#round depth to nearest integer
start_dates$depth <- ceiling(start_dates$depth)
params$depth <- ceiling(params$depth)

#FIGS
ggplot(subset(start_dates, depth %in% c(1,5,9)),
       aes(DA, RMSE, fill=as.factor(da_start_date))) +  ylab("RMSE") + xlab("")+
  geom_bar(stat="identity",position="dodge") + theme_bw() + guides(fill=guide_legend(title="")) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = c(0.75,0.31),
        legend.background = element_blank(),legend.direction = "horizontal", panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0.05,-0.2,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 6),legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  facet_grid(depth~phen, scales="free",labeller = labeller(depth = depths)) + scale_fill_manual(values=c("#81A665","#E0CB48","#D08151")) 
ggsave(file.path(lake_directory,"analysis/figures/RMSEvsDAfreq_depth_facets_start_dates.jpg"),width=3.5, height=4)

ggplot(subset(params, depth %in% c(1,5,9)),
       aes(DA, RMSE, fill=as.factor(param))) +  ylab("RMSE") + xlab("")+
  geom_bar(stat="identity",position="dodge") + theme_bw() + guides(fill=guide_legend(title="Params")) +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = c(0.75,0.31),
        legend.background = element_blank(),legend.direction = "horizontal", panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0.05,-0.2,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 6),legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  facet_grid(depth~phen, scales="free",labeller = labeller(depth = depths)) + scale_fill_manual(values=c("#81A665","#E0CB48","#D08151")) 
ggsave(file.path(lake_directory,"analysis/figures/RMSEvsDAfreq_depth_facets_parameters.jpg"),width=3.5, height=4)

#------------------------------------------------------------------------------------------------#
#parameter evolution figs 
#read in all forecasts 
params_dir <- arrow::SubTreeFileSystem$create(file.path(lake_directory,"scores/da_study"))
params <- arrow::open_dataset(score_dir) |> collect() |>   
  filter(variable %in% c("lw_factor","zone1temp","zone2temp"), horizon >=0)

#change datetime format
params$datetime <- as.Date(params$datetime)

#change model_id to be all uppercase
params$model_id <- str_to_title(params$model_id)

#change DA factor order
params$model_id <- factor(params$model_id, levels = c("Daily", "Weekly","Fortnightly","Monthly"))

ggplot(subset(params,horizon==0), aes(datetime, mean, color=model_id)) + theme_bw() +
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
#ggsave(file.path(lake_directory,"analysis/figures/paramevolvstime_0day.jpg"),width=3.5, height=4)

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



source(file.path(lake_directory,"R/old/read_flare_params.R"))

forecasts_daily_nc <- list.files(file.path(lake_directory,"forecasts/bvre/DA_experiments/27_nov_start/daily"), pattern=".nc", full.names=TRUE)[-c(1)] #ignoring first file because this is the DA period
forecasts_weekly_nc <- list.files(file.path(lake_directory,"forecasts/bvre/DA_experiments/27_nov_start/weekly"), pattern=".nc", full.names=TRUE)[-c(1)]
forecasts_fortnightly_nc <- list.files(file.path(lake_directory,"forecasts/bvre/DA_experiments/27_nov_start/fortnightly"), pattern=".nc", full.names=TRUE)[-c(1)]
forecasts_monthly_nc <- list.files(file.path(lake_directory,"forecasts/bvre/DA_experiments/27_nov_start/monthly"), pattern=".nc", full.names=TRUE)[-c(1)]

forecasts_daily_nc_24nov <- list.files(file.path(lake_directory,"forecasts/bvre/DA_experiments/24_nov_start/daily"), pattern=".nc", full.names=TRUE)[-c(1)] #ignoring first file because this is the DA period
forecasts_weekly_nc_24nov <- list.files(file.path(lake_directory,"forecasts/bvre/DA_experiments/24_nov_start/weekly"), pattern=".nc", full.names=TRUE)[-c(1)]
forecasts_fortnightly_nc_24nov <- list.files(file.path(lake_directory,"forecasts/bvre/DA_experiments/24_nov_start/fortnightly"), pattern=".nc", full.names=TRUE)[-c(1)]
forecasts_monthly_nc_24nov <- list.files(file.path(lake_directory,"forecasts/bvre/DA_experiments/24_nov_start/monthly"), pattern=".nc", full.names=TRUE)[-c(1)]

forecasts_daily_nc_22nov <- list.files(file.path(lake_directory,"forecasts/bvre/DA_experiments/22_nov_start/daily"), pattern=".nc", full.names=TRUE)[-c(1)] #ignoring first file because this is the DA period
forecasts_weekly_nc_22nov <- list.files(file.path(lake_directory,"forecasts/bvre/DA_experiments/22_nov_start/weekly"), pattern=".nc", full.names=TRUE)[-c(1)]
forecasts_fortnightly_nc_22nov <- list.files(file.path(lake_directory,"forecasts/bvre/DA_experiments/22_nov_start/fortnightly"), pattern=".nc", full.names=TRUE)[-c(1)]
forecasts_monthly_nc_22nov <- list.files(file.path(lake_directory,"forecasts/bvre/DA_experiments/22_nov_start/monthly"), pattern=".nc", full.names=TRUE)[-c(1)]

forecasts_daily_nc_1.04 <- list.files(file.path(lake_directory,"forecasts/bvre/DA_experiments/inflat_1.04/daily"), pattern=".nc", full.names=TRUE)[-c(1)] #ignoring first file because this is the DA period
forecasts_weekly_nc_1.04 <- list.files(file.path(lake_directory,"forecasts/bvre/DA_experiments/inflat_1.04/weekly"), pattern=".nc", full.names=TRUE)[-c(1)]
forecasts_fortnightly_nc_1.04 <- list.files(file.path(lake_directory,"forecasts/bvre/DA_experiments/inflat_1.04/fortnightly"), pattern=".nc", full.names=TRUE)[-c(1)]
forecasts_monthly_nc_1.04 <- list.files(file.path(lake_directory,"forecasts/bvre/DA_experiments/inflat_1.04/monthly"), pattern=".nc", full.names=TRUE)[-c(1)]

forecasts_daily_nc_fixed <- list.files(file.path(lake_directory,"forecasts/bvre/DA_experiments/fixed_params/daily"), pattern=".nc", full.names=TRUE)[-c(1)] #ignoring first file because this is the DA period
forecasts_weekly_nc_fixed <- list.files(file.path(lake_directory,"forecasts/bvre/DA_experiments/fixed_params/weekly"), pattern=".nc", full.names=TRUE)[-c(1)]
forecasts_fortnightly_nc_fixed <- list.files(file.path(lake_directory,"forecasts/bvre/DA_experiments/fixed_params/fortnightly"), pattern=".nc", full.names=TRUE)[-c(1)]
forecasts_monthly_nc_fixed <- list.files(file.path(lake_directory,"forecasts/bvre/DA_experiments/fixed_params/monthly"), pattern=".nc", full.names=TRUE)[-c(1)]



#summary stats for first forecast of each
daily <- read_flare_params(files = forecasts_daily_nc, type = "forecast", summary = TRUE) %>% mutate(DA="Daily")
weekly <- read_flare_params(files = forecasts_weekly_nc, type = "forecast", summary = TRUE) %>% mutate(DA="Weekly")
fortnightly <- read_flare_params(files = forecasts_fortnightly_nc, type = "forecast", summary = TRUE) %>% mutate(DA="Fortnightly")
monthly <- read_flare_params(files = forecasts_monthly_nc, type = "forecast", summary = TRUE) %>% mutate(DA="Monthly")

daily_24nov <- read_flare_params(files = forecasts_daily_nc_24nov, type = "forecast", summary = TRUE) %>% mutate(DA="Daily")
weekly_24nov <- read_flare_params(files = forecasts_weekly_nc_24nov, type = "forecast", summary = TRUE) %>% mutate(DA="Weekly")
fortnightly_24nov <- read_flare_params(files = forecasts_fortnightly_nc_24nov, type = "forecast", summary = TRUE) %>% mutate(DA="Fortnightly")
monthly_24nov <- read_flare_params(files = forecasts_monthly_nc_24nov, type = "forecast", summary = TRUE) %>% mutate(DA="Monthly")

daily_22nov <- read_flare_params(files = forecasts_daily_nc_22nov, type = "forecast", summary = TRUE) %>% mutate(DA="Daily")
weekly_22nov <- read_flare_params(files = forecasts_weekly_nc_22nov, type = "forecast", summary = TRUE) %>% mutate(DA="Weekly")
fortnightly_22nov <- read_flare_params(files = forecasts_fortnightly_nc_22nov, type = "forecast", summary = TRUE) %>% mutate(DA="Fortnightly")
monthly_22nov <- read_flare_params(files = forecasts_monthly_nc_22nov, type = "forecast", summary = TRUE) %>% mutate(DA="Monthly")

daily_1.04 <- read_flare_params(files = forecasts_daily_nc_1.04, type = "forecast", summary = TRUE) %>% mutate(DA="Daily")
weekly_1.04 <- read_flare_params(files = forecasts_weekly_nc_1.04, type = "forecast", summary = TRUE) %>% mutate(DA="Weekly")
fortnightly_1.04 <- read_flare_params(files = forecasts_fortnightly_nc_1.04, type = "forecast", summary = TRUE) %>% mutate(DA="Fortnightly")
monthly_1.04 <- read_flare_params(files = forecasts_monthly_nc_1.04, type = "forecast", summary = TRUE) %>% mutate(DA="Monthly")

#fixed parameters --> zone1temp set to 5.54


#combine all parameter dfs
parameters <- rbind(daily, weekly, fortnightly, monthly)
parameters$inflat <- 1.02

parameters_1.04 <- rbind(daily_1.04, weekly_1.04, fortnightly_1.04, monthly_1.04)
parameters_1.04$inflat <- 1.04

params_inflat <- rbind(parameters,parameters_1.04)

parameters_24nov <- rbind(daily_24nov, weekly_24nov, fortnightly_24nov, monthly_24nov)
parameters_22nov <- rbind(daily_22nov, weekly_22nov, fortnightly_22nov, monthly_22nov)

parameters_24nov$inflat <- 1.02
parameters_22nov$inflat <- 1.02

parameters$start_date <- "27nov"
parameters_24nov$start_date <- "24nov"
parameters_22nov$start_date <- "22nov"

params_startdate <- rbind(parameters,parameters_22nov,parameters_24nov)

#change DA factor order
params_inflat$DA <- factor(params_inflat$DA, levels = c("Daily", "Weekly","Fortnightly","Monthly"))
params_startdate$DA <- factor(params_startdate$DA, levels = c("Daily", "Weekly","Fortnightly","Monthly"))

#rename zone1temp so more informative facet label in fig
params_inflat <- params_inflat %>% mutate(parameter = recode(parameter, "zone1temp" = "Sediment Temperature"))
params_startdate <- params_startdate %>% mutate(parameter = recode(parameter, "zone1temp" = "Sediment Temperature"))

#visualize how parameters change over time
ggplot(subset(params_inflat, parameter=="Sediment Temperature"), aes(datetime, mean, color=DA)) + theme_bw() +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = c(0.59,0.03),
        legend.background = element_blank(),legend.direction = "horizontal", panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0.05,-0.2,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 6),legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  scale_color_manual(values=cb_friendly_2) +
  facet_wrap(~inflat, scales="free_y", ncol=1)  + scale_fill_manual(values=cb_friendly_2) +
  scale_x_date(date_labels = "%b") + ylab(expression("Temperature ("*~degree*C*")")) + xlab("")  +
  geom_ribbon(aes(y = mean, ymin = mean-sd, ymax = mean+sd, color=DA, fill=DA), alpha=0.5) +
  guides(fill = guide_legend(title="DA frequency"), color = guide_legend(title="DA frequency"))
ggsave(file.path(lake_directory,"analysis/figures/paramRMSEvsHorizon_inflat.jpg"),width=3.5, height=4)

ggplot(subset(params_startdate, parameter=="Sediment Temperature"), aes(datetime, mean, color=DA)) + theme_bw() +
  theme(text = element_text(size=8), axis.text = element_text(size=6, color="black"), legend.position = c(0.59,0.03),
        legend.background = element_blank(),legend.direction = "horizontal", panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0.05,-0.2,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 6),legend.text  = element_text(size = 6), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6), axis.text.y = element_text(size=6)) +
  scale_color_manual(values=cb_friendly_2) +
  facet_wrap(~start_date, scales="free_y", ncol=1)  + scale_fill_manual(values=cb_friendly_2) +
  scale_x_date(date_labels = "%b") + ylab(expression("Temperature ("*~degree*C*")")) + xlab("")  +
  geom_ribbon(aes(y = mean, ymin = mean-sd, ymax = mean+sd, color=DA, fill=DA), alpha=0.5) +
  guides(fill = guide_legend(title="DA frequency"), color = guide_legend(title="DA frequency"))
ggsave(file.path(lake_directory,"analysis/figures/paramRMSEvsHorizon_start_date.jpg"),width=3.5, height=4)

#--------------------------------------------------------------------------------------------#
# Figure comparing Mixed w/ ice-cover data and Mixed w/o ice-cover data
# ice-on/off dates for BVR 2021: 10Jan/12Jan, 30Jan/31Jan 13Feb/16Feb

#creating smaller dataset for kw test w/ 1,5,9m and 1,7,35 days
kw_horizons <- forecast_skill_depth_horizon_27nov[forecast_skill_depth_horizon_27nov$depth %in% c(1,5,9) & forecast_skill_depth_horizon_27nov$horizon %in% c(1,7,35) & 
                                                    forecast_skill_depth_horizon_27nov$DA %in% c("Daily","Weekly","Fortnightly","Monthly"),]

#only select mixed period
kw_horizons_mixed <- kw_horizons[kw_horizons$phen=="Mixed",]


#create new df with all mixed days AND mixed days w/o ice
kw_horizons_mixed_sub <-   rbind(
  cbind(kw_horizons_mixed, faceter = "all"),
  cbind(kw_horizons_mixed[!(kw_horizons_mixed$forecast_date %in% 
                              c(as.Date("2021-01-10"), as.Date("2021-01-11"),as.Date("2021-01-30"),
                                as.Date("2021-02-13"),as.Date("2021-02-14"),as.Date("2021-02-15"))),],
        faceter = "no ice")
)


#rename depth and ice facets
faceter <- c("Mixed with ice","Mixed without ice")
names(faceter) <- c("all","no ice")

depths <- c("1m","5m","9m")
names(depths) <- c("1","5","9")

#order factor levels
kw_horizons_mixed_sub$DA <- factor(kw_horizons_mixed_sub$DA, levels = c("Daily", "Weekly", "Fortnightly", "Monthly"))

kw_horizons_mixed_sub %>%
  group_by(DA,depth,horizon,faceter) %>%  # do the same calcs for each box
  mutate(value2 = filter_lims(RMSE)) %>%
  ggplot(aes(DA, value2, fill=as.factor(horizon))) +  ylab("RMSE") + xlab("")+
  geom_boxplot(outlier.shape = NA) + theme_bw() + guides(fill=guide_legend(title="Horizon (days)")) +
  geom_hline(yintercept=2, linetype='dashed', col = 'black') +
  theme(text = element_text(size=4), axis.text = element_text(size=6, color="black"), legend.position = c(0.77,0.24),
        legend.background = element_blank(),legend.direction = "horizontal", panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0.05,-0.2,0), "cm"),legend.key.size = unit(0.5, "lines"), panel.grid.major = element_blank(),
        legend.title = element_text(size = 3),legend.text  = element_text(size = 3), panel.spacing=unit(0, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=4), axis.text.y = element_text(size=4)) +
  #geom_text(data=letters,aes(x=DA,y=0.2+max.RMSE,label=letters$letter),hjust=0.1,vjust = -0.1, size=1.5) +
  facet_grid(depth~faceter, scales="free_y",labeller = labeller(faceter = faceter, depth = depths)) + scale_fill_manual(values=c("#81A665","#E0CB48","#D08151")) 
ggsave(file.path(lake_directory,"analysis/figures/RMSEvsDAfreq_depth_facets_IcevsNoice.jpg"))



#2021 phenology: 2021-03-08 is first time when >3 consecutive days had difference between surface and bottom >1C
#code for calculating strat/mixed periods
# bvr_temps <- temp_long %>% filter(temp_long$Variable=="temperature" & DateTime>= "2021-01-01") %>% select(DateTime, Reading, Depth) %>%
#   mutate(DateTime = as.Date(DateTime))  %>% mutate(Depth = round(Depth,0))
# 
# 
# bvr_surf_bot_temps <- bvr_temps %>% group_by(DateTime, Depth) %>% summarise(Temp = mean(Reading)) %>% 
#   group_by(DateTime)  %>% filter(Depth== min(Depth) | Depth== max(Depth)) %>%
#   mutate(diff = abs(last(Temp)-first(Temp))) %>% mutate(diff = round(diff,0))
# 
# plot(bvr_surf_bot_temps$DateTime, bvr_surf_bot_temps$diff)
# 
# mix <- bvr_surf_bot_temps[bvr_surf_bot_temps$diff<1,] 
# strat <- bvr_surf_bot_temps[bvr_surf_bot_temps$diff>=1,]

