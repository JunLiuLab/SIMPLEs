library(matrixStats)
library(here)
library(tidyverse)


## ref: github repo: imputationBenchmark
## data/code/process/cellbench/04_5celline.R

## * load data
## data from github repo: sc_mixology
## data/csv/
## need to gunzip the csv.zip files.

cnt <- read.csv(
  file = here::here("10xGenomics", "scRNAseq", "csv", "sc_10x_5cl.count.csv"),
  header = TRUE, row.names = 1
)

scmeta <- read.csv(
  file = here::here(
    "10xGenomics", "scRNAseq", "csv",
    "sc_10x_5cl.metadata.csv"
  ),
  header = TRUE, row.names = 1
)
## match column with cell line names.
colnames(cnt) <- paste0(
  colnames(cnt), ":",
  scmeta[match(colnames(cnt), rownames(scmeta)), "cell_line_demuxlet"]
)

## * clean cnt data.

## ** keep genes expressed at least in 10% cells
## nec short for number of expressed cells for gene
necg <- rowSums(cnt > 0)
tc <- dim(cnt)[2]
rgcnt <- cnt[which(necg > 0.1 * tc), ]

## ** remove mitochondiron genes
## mm short for move mitochondiron
rgmmcnt <- rgcnt[!grepl("^MT-", rownames(rgcnt)), ]

## * save data
saveRDS(rgmmcnt, file = here::here("10xGenomics", "scRNAseq", "sc_10x_5cl_rgmm_cnt.rds"))
write.csv(rgmmcnt,
  file = here::here("10xGenomics", "scRNAseq", "sc_10x_5cl_rgmm_cnt.csv"),
  quote = FALSE, row.names = TRUE
)
