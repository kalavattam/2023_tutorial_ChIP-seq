#!/usr/bin/env Rscript

#  set-up_conda-R_2023-0518.R
#  KA

#  On setting up an R environment via conda/mamba
#  1. Set up an environment for doing work with R packages; on the command line:
#+    mamba create -n R_env -c conda-forge r-base==4.2.3 r-biocmanager==1.30.20 r-tidyverse
#+    
#+ 2. Open RStudio in such a way that it's "within" the above environment, R_env; on the command line:
#+    mamba activate R_env && open -na /Applications/RStudio.app  # It's helpful to make an alias for this
#+    
#+ 3. Install Bioconductor packages not through conda/mamba but instead through BiocManager, a program for managing Bioconductor packages via R

suppressMessages(library(BiocManager))

if(!require("csaw", quietly = TRUE))
    BiocManager::install("csaw")
