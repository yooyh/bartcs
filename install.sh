#!/bin/bash

set -e

## build ARGs
NCPUS=${NCPUS:--1}

# install packages for bartcs
install2.r --error --skipinstalled -n "$NCPUS" \
    devtools \
    ggcharts \
    ggplot2 \
    invgamma \
    MCMCpack \
    Rcpp \
    rlang \
    rootSolve \
    stats \
    knitr \
    microbenchmark \
    rmarkdown \
    attachment \
    checkhelper \
    rhub \
    spelling

# Clean up
rm -rf /var/lib/apt/lists/*
rm -rf /tmp/downloaded_packages