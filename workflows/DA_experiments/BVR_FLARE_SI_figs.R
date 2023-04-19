#BVR FLARE SI figs

#load libraries
pacman::p_load(reshape2, egg, dplyr, stringr)

#change tag_facet code
tag_facet2 <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                       hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
}

#------------------------------------------------------------------------------------------------#
#parameter evolution figs 
#read in all forecasts 
params_dir <- arrow::SubTreeFileSystem$create(file.path(lake_directory,"scores/all_UC"))
params <- arrow::open_dataset(params_dir) |> collect() |>   
  filter(variable %in% c("lw_factor","zone1temp","zone2temp"), horizon >=0)

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
ggsave(file.path(lake_directory,"analysis/figures/paramevolvstime_1day_fig5.jpg"),width=3.5, height=4)

#figuring out the date that DA parameters diverge
mean(params$mean[params$variable=="lw_factor" & params$model_id=="Daily" & params$datetime >= "2021-04-01"])

mean(params$mean[params$variable=="lw_factor" & params$model_id=="Weekly" & params$datetime >= "2021-04-01"])
mean(params$mean[params$variable=="lw_factor" & params$model_id=="Fortnightly" & params$datetime >= "2021-04-01"])
mean(params$mean[params$variable=="lw_factor" & params$model_id=="Monthly" & params$datetime >= "2021-04-01"])

range(c(mean(params$mean[params$variable=="zone2temp" & params$model_id=="Weekly" & params$datetime >= "2021-04-01"]),
       mean(params$mean[params$variable=="zone2temp" & params$model_id=="Fortnightly" & params$datetime >= "2021-04-01"]),
       mean(params$mean[params$variable=="zone2temp" & params$model_id=="Monthly" & params$datetime >= "2021-04-01"])))

range(params$mean[params$variable=="zone2temp" & params$model_id!="Daily" & params$datetime >= "2021-04-01"])
range(params$mean[params$variable=="zone1temp" & params$model_id!="Daily" & params$datetime >= "2021-04-01"])

range(c(mean(params$mean[params$variable=="zone1temp" & params$model_id=="Weekly" & params$datetime >= "2021-04-01"]),
        mean(params$mean[params$variable=="zone1temp" & params$model_id=="Fortnightly" & params$datetime >= "2021-04-01"]),
        mean(params$mean[params$variable=="zone1temp" & params$model_id=="Monthly" & params$datetime >= "2021-04-01"])))

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

