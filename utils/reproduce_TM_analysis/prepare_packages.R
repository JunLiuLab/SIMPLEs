install.packages("devtools")
devtools::install_github("klutometis/roxygen")

# install SIMPLEs
devtools::install_github("JunLiuLab/SIMPLEs", ref = "master", build_vignettes = TRUE)

install.packages("Seurat")

# For directory handle
install.packages("here")
devtools::install_github("smbache/import")

# grammar of data manipulation on data.frame
install.packages("tidyverse")

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("DESeq2")

install.packages("Matrix")
install.packages("matrixStats")
install.packages("ontologyIndex")
install.packages("ggplot2")
install.packages("doParallel")
install.packages("optparse")

# for drawing images
devtools::install_github("kassambara/ggpubr")
install.packages("gridExtra")
install.packages("gridGraphics")
install.packages("cowplot")
