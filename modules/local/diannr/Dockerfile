# FROM continuumio/miniconda3:latest
FROM condaforge/mambaforge:latest
# USER root

RUN mamba install --yes -c conda-forge r-base r-essentials r-devtools r-argparse

RUN R -e "devtools::install_github('https://github.com/vdemichev/diann-rpackage')"
