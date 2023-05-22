FROM rocker/geospatial:4.1.2

RUN apt-get update && apt-get -y install libgd-dev libnetcdf-dev git

RUN git clone https://github.com/EcoDynForecast/BVRE-forecast-code.git /home/rstudio/BVRE-forecast-code



