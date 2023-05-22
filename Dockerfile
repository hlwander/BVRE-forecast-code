FROM rocker/geospatial:4.1.2

RUN apt-get update && apt-get -y install libgd-dev libnetcdf-dev git

USER rstudio

RUN git clone https://github.com/EcoDynForecast/BVRE-forecast-code.git /home/rstudio/BVRE-forecast-code

#RUN chmod 777 /home/rstudio/BVRE-forecast-code

#RUN Rscript /home/rstudio/BVRE-forecast-code/workflows/DA_experiments/install.R


