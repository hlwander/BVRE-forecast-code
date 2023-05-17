#load packages
if (!require("pacman"))install.packages("pacman")
pacman::p_load(curl, raster, lubridate, tidyverse, zoo, aws.s3) 

# install packages
install.packages(c("remotes", "tidyverse"))
install.packages(c("here", "aws.s3"))
remotes::install_github("SwampThingPaul/AnalystHelper")
remotes::install_github("FLARE-forecast/GLM3r") 
remotes::install_github("rqthomas/glmtools", ref = "b50e9a7b73e41afcd8119e2b9ac172c1c7beb51f")
remotes::install_github("rqthomas/FLAREr", ref = "40b9bfd6a4ec54942364cf6eed2aa46f022eeefa")
