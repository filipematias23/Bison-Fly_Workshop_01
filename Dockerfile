FROM rocker/binder:4.3.1

# Copy repo into ${HOME}, make user own $HOME
USER root
COPY . ${HOME}
RUN chown -R ${NB_USER} ${HOME}
USER ${NB_USER}
  
## run any install.R script we find
RUN if [ -f install.R ]; then R --quiet -f install.R; fi
