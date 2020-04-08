## Base image
FROM r-base

# Install R devtools
RUN apt-get update && \
  apt-get install -y libcurl4-openssl-dev libssl-dev libssh2-1-dev libxml2-dev && \
  R -e "install.packages('remotes', dependencies = c('Imports', 'Depends'))"

# Install mHMMbayes
RUN R -e "remotes::install_github('emmekeaarts/mHMMbayes@369420e05cb96e6af761a8bd9bd30b8539d27a24', dependencies='Imports')"
RUN R -e "install.packages('argparser')"
RUN R -e "install.packages('logger')"
RUN R -e "install.packages('uuid')"
RUN R -e "remotes::install_github('JasperHG90/sleepsimR@99e7a8e0dabf4f732d0ce1bab3a8a527f42ca73d', dependencies='Imports')"

## Copy app
COPY ./app /app

## Set entry point
ENTRYPOINT ["Rscript", "app/run_sleepdata_analysis.R"]
