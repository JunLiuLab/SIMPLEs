library(doParallel)
library(optparse)
library(matrixStats)
library(here)
library(tidyverse)
## * methods to compare
library(SIMPLEs)
library(SCRABBLE)
library(SAVER)
library(Rmagic)
library(VIPER)
library(scImpute)


## * load option
option_list <- list(
  make_option(c("-m", "--method"),
    action = "store",
    type = "character", default = "control"
  ),
  make_option(c("--cnm"),
    action = "store",
    type = "integer", default = 1000
    ),
  make_option(c("-n", "--ncores"),
              action = "store",
              type = "integer", default = 1
              )
)
args <- option_list %>%
  OptionParser(option_list = .) %>%
  parse_args()

get_tpm <- function(obs_counts_gbc, scale_level=1.0) {
  scale(obs_counts_gbc, center = F, scale = colSums(obs_counts_gbc)) *
    scale_level
}

log_tpm_rseq <- function(tpm_gbc, margin = 1) {
  tmp <- log(tpm_gbc + margin)
  ## add an small value to first posion if that row are equals
  ## to the same value.
  eps <- 0.1
  tmp_rowvar <- matrixStats::rowVars(tmp)
  modify_rows <- which(dplyr::near(tmp_rowvar, 0))
  tmp[modify_rows, 1] <- tmp[modify_rows, 1] + eps
  return(tmp)
}


load_sampled_data <- function(mycnm) {
  fnm <- paste0("leu_simple_cell_", mycnm, ".rds")
  readRDS(here::here("TabulaMuris", "supplement", "timer", fnm))
}


## * impute with different methods.
control_impute <- function(logtpm, ...) {
  return(logtpm)
}

simple_impute <- function(logtpm, k = 10, m = 1, mingene = 100,
                          max_lambda = T, pm = 0.4, cutoff = 0.1, ...) {
  r <- SIMPLEs::SIMPLE(
    dat = logtpm, K0 = k, M0 = m, max_lambda = max_lambda,
    min_gene = mingene, p_min = pm, cutoff = cutoff
  )
  return(r$impt)
}

saver_impute <- function(tpm, ncores = 1, ...) {
  r <- SAVER::saver(x = tpm, ncores = ncores, size.factor = 1)
  return(log(r$estimate + 1))
}

## must library(VIPER)
viper_impute <- function(tpm, ...) {
  num <- as.integer(dim(tpm)[1] * 0.7)
  r <- VIPER::VIPER(
    gene.expression = tpm, num = num, percentage.cutoff = 0.1,
    minbool = FALSE, alpha = 0.5, report = F, outdir = NULL,
    prefix = NULL
  )
  return(r$imputed_log)
}

scrabble_impute <- function(logtpm, ...) {
  SCRABBLE::scrabble(
    data = list(data_sc = logtpm, data_bulk = NULL),
    parameter = c(1, 1e-06, 1e-04), nIter = 20,
    error_out_threshold = 1e-04,
    nIter_inner = 20, error_inner_threshold = 1e-04
  )
}

scimpute_impute <- function(obs_count, scale_level, ncores = 1, ...) {
  remove_temp_file <- function() {
    if (file.exists(mytempfilenm)) {
      file.remove(mytempfilenm)
    }
    if (file.exists(rsf)) {
      file.remove(rsf)
    }
  }
  rownames(obs_count) <- seq_len(nrow(obs_count))
  colnames(obs_count) <- seq_len(ncol(obs_count))

  mytempfilenm <- paste0("scImpute", sample.int(10000, 1), sample.int(1000, 1), "input", ".csv")
  write.csv(obs_count, file = mytempfilenm, row.names = T)

  scImpute::scimpute(
    count_path = mytempfilenm, infile = "csv",
    outfile = "csv", out_dir = "./", Kcluster = 2, ncores = ncores,
    drop_thre = 0.5
  )
  rsf <- "scimpute_count.csv"
  if (file.exists(rsf)) {
    imp_count <- read.csv(file = rsf, header = TRUE, row.names = 1)
    tpm <- get_tpm(imp_count, scale_level)
    remove_temp_file()
    return(log(tpm + 1))
  } else {
    remove_temp_file()
    stop("scImpute did not output result.")
  }
}

rmagic_impute <- function(logtpm, ...) {
  ## rmagic nee cell by gene matrix, while logtpm is gene by
  ## cell
  imp_cbg <- magic(t(logtpm),
    genes = "all_genes", knn = 5,
    t = "auto"
  )
  return(t(imp_cbg$result))
}

## * run method and evaluate it.
main_impute <- function() {
  logtpm_tm <- load_sampled_data(args$cnm)
  tpm_tm <- exp(logtpm_tm) - 1.0
  near_cnt <- round(tpm_tm)
  k <- 10
  m <- 1

  ## ** parameters
  method <- args$method
  ## ** imputation
  mytimer <- system.time({
    do.call(
      what = str_c(args$method, "_impute"),
      args = list(
        method = method, tpm = tpm_tm, logtpm = logtpm_tm,
        k = k, m = m, ncores = args$ncores, obs_count = near_cnt,
        scale_level = 1.0
      )
    )
  })
  print(mytimer)


  saveRDS(object = mytimer, file = here::here(
    "TabulaMuris", "supplement", "timer",
    paste0(method, "_timer_", args$cnm, ".rds")
  ))
}
## * script main interface
registerDoParallel(cores = args$ncores)
main_impute()
