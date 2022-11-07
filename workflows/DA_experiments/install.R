# install packages
install.packages(c("remotes", "tidyverse"))
install.packages(c("here", "aws.s3"))
remotes::install_github("FLARE-forecast/GLM3r", ref = "6561ffd5ee5807dd1cac5cc9bb0fff7fc4c43b0b")
remotes::install_github("rqthomas/glmtools", ref = "b50e9a7b73e41afcd8119e2b9ac172c1c7beb51f")
remotes::install_github("rqthomas/FLAREr", ref = "1b705ca7bd19f1fdc0d8befe21a91eee6d19d3c0")
