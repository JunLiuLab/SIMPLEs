install.packages("devtools")
library(devtools)
devtools::install_github("klutometis/roxygen")
library(roxygen2)

devtools::install_github("YosefLab/SymSim")

devtools::install_github("JunLiuLab/SIMPLEs",
  ref = "master",
  build_vignettes = TRUE
)

devtools::install_github("Vivianstats/scImpute")
devtools::install_github("mohuangx/SAVER")
devtools::install_github("software-github/SCRABBLE/R")

install.packages("Rcpp")
install.packages("RcppArmadillo")
install.packages("quadprog")
install.packages("glmnet")
devtools::install_github("ChenMengjie/VIPER")

## NOTE: pip install --user magic-impute <---- is needed.
install.packages("Rmagic")

# For directory handle
install.packages("here")
library(here)

# grammar of data manipulation on data.frame
install.packages("tidyverse")
library(dplyr)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install("DESeq2")

devtools::install_github("smbache/import")

install.packages("Matrix")
install.packages("matrixStats")
install.packages("ontologyIndex")
install.packages("ggplot2")
install.packages("doParallel")
install.packages("optparse")
## or install optparse below.
devtools::install_github("trevorld/r-getopt")
devtools::install_github("trevorld/r-optparse")
