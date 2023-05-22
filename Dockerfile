FROM rocker/geospatial:4.1.2

## Declares build arguments
ARG NB_USER
ARG NB_UID

COPY --chown=${NB_USER} . ${HOME}

ENV DEBIAN_FRONTEND=noninteractive
USER root
RUN echo "Checking for 'apt.txt'..." \
        ; if test -f "apt.txt" ; then \
        apt-get update --fix-missing > /dev/null\
        && xargs -a apt.txt apt-get install --yes \
        && apt-get clean > /dev/null \
        && rm -rf /var/lib/apt/lists/* \
        ; fi
USER ${NB_USER}

## Run an install.R script, if it exists.
#RUN if [ -f workflows/DA_experiments/install.R ]; then R --quiet -f workflows/DA_experiments/install.R; fi

RUN mkdir /home/rstudio/wander_et_al      
RUN cd /home/rstudio/wander_et_al
RUN pwd
RUN git clone https://github.com/EcoDynForecast/BVRE-forecast-code.git

