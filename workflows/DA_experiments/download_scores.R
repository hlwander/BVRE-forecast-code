lake_directory <- here::here()
download.file(url = "https://zenodo.org/record/7951402/files/scores.zip?download=1", destfile = file.path(lake_directory,"scores.zip"), method = "curl")
unzip(file.path(lake_directory,"scores.zip"))

download.file(url = "https://zenodo.org/record/7951402/files/targets.zip?download=1", destfile = file.path(lake_directory,"targets.zip"), method = "curl")
unzip(file.path(lake_directory,"targets.zip"))
