# FROM continuumio/miniconda3:latest
FROM condaforge/mambaforge:latest
# USER root

RUN mamba install --yes -c conda-forge r-base r-essentials r-devtools r-argparse

RUN Rscript -e "install.packages('iq',dependencies=TRUE, repos='http://cran.rstudio.com/')"
