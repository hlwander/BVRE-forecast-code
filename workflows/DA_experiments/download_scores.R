lake_directory <- here::here()
download.file(url = "https://zenodo.org/record/7925098/files/scores.zip?download=1", destfile = file.path(lake_directory,"archive/scores.zip"), method = "curl")
unzip(file.path(lake_directory,"archive/scores.zip") ,exdir = "scores")

download.file(url = "https://zenodo.org/record/7925098/files/targets.zip?download=1", destfile = file.path(lake_directory,"archive/targets.zip"), method = "curl")
unzip(file.path(lake_directory,"archive/targets.zip") ,exdir = "targets")
