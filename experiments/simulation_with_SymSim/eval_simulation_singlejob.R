library(doParallel)
library(optparse)
library(matrixStats)
library(caTools)
library(here)
library(tidyverse)
## * methods to compare
library(SIMPLEs)
library(SCRABBLE)
library(VIPER)
library(SAVER)
library(Rmagic)
library(scImpute)

## * get optparse ref do_imputation in TM analysis.
option_list <- list(
  make_option(c("-m", "--method"),
    action = "store",
    type = "character", default = "control"
  ),
  make_option(c("-k", "--K0"),
    action = "store", type = "integer",
    default = 2
  ),
  make_option(c("--M0"),
    action = "store", type = "integer",
    default = 1
  ),
  make_option(c("--nrep"),
    action = "store", type = "integer",
    default = 10
  )
)
args <- option_list %>%
  OptionParser(option_list = .) %>%
  parse_args()
## * load data
## detest could be t.test or wilcox test
p <- list(
  nrep = args$nrep, npop = 5, ncell = 300, ngene = 1000, ncores = 8, scale_level_umi = 10000,
  scale_level_nonumi = 1e+06, method = args$method, k = args$K0,
  m = args$M0, pcs = seq(2, 10, 2), detest = "t"
)

## gbc short for gene by cell
get_tpm <- function(obs_counts_gbc, scale_level) {
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

symsim_get_degenes <- function(true_counts_res, popA, popB) {
  meta_cell <- true_counts_res$cell_meta
  meta_gene <- true_counts_res$gene_effects
  popA_idx <- which(meta_cell$pop == popA)
  popB_idx <- which(meta_cell$pop == popB)
  ngenes <- dim(true_counts_res$gene_effects[[1]])[1]

  DEstr <- sapply(strsplit(colnames(meta_cell)[which(grepl("evf", colnames(meta_cell)))], "_"), "[[", 2)
  param_str <- sapply(strsplit(colnames(meta_cell)[which(grepl("evf", colnames(meta_cell)))], "_"), "[[", 1)
  n_useDEevf <- sapply(1:ngenes, function(igene) {
    return(sum(abs(meta_gene[[1]][igene, DEstr[which(param_str == "kon")] == "DE"]) - 0.001 > 0) +
      sum(abs(meta_gene[[2]][igene, DEstr[which(param_str == "koff")] == "DE"]) - 0.001 > 0) +
      sum(abs(meta_gene[[3]][igene, DEstr[which(param_str == "s")] == "DE"]) - 0.001 > 0))
  })

  kon_mat <- true_counts_res$kinetic_params[[1]]
  koff_mat <- true_counts_res$kinetic_params[[2]]
  s_mat <- true_counts_res$kinetic_params[[3]]

  logFC_theoretical <- sapply(1:ngenes, function(igene) {
    return(log2(mean(s_mat[igene, popA_idx] * kon_mat[igene, popA_idx] / (kon_mat[igene, popA_idx] + koff_mat[igene, popA_idx])) /
      mean(s_mat[igene, popB_idx] * kon_mat[igene, popB_idx] / (kon_mat[igene, popB_idx] + koff_mat[igene, popB_idx]))))
  })

  true_counts <- true_counts_res$counts
  true_counts_norm <- t(t(true_counts) / colSums(true_counts)) * 10^6

  p_true_counts <- sapply(1:ngenes, function(igene) {
    return(wilcox.test(true_counts_norm[igene, popA_idx], true_counts_norm[igene, popB_idx])$p.value)
  })

  adjp_true_counts <- p.adjust(p_true_counts, method = "fdr")

  return(list(nDiffEVF = n_useDEevf, logFC_theoretical = logFC_theoretical, wil.p_true_counts = adjp_true_counts))
}


load_symsim_data <- function(myseed) {
  ## load obs_nonumi, obs_umi, and true_counts
  load(file = str_c("symsim_data/sim", p$ncell, p$ngene, myseed, ".RData", sep = "_"))
  tpm_gbc <- get_tpm(obs_nonumi$counts, p$scale_level_nonumi)
  logtpm_gbc <- log_tpm_rseq(tpm_gbc)
  ## tth is short for ten thousand
  tptth_gbc <- get_tpm(obs_umi$counts, p$scale_level_umi)
  logtptth_gbc <- log_tpm_rseq(tptth_gbc)
  ## ** clusters cc short for cell cluster
  cc <- true_counts$cell_meta$pop
  npop <- nlevels(as.factor(cc))
  ## ** get ground truths of DE(differential expression) genes
  ## and cell clusters.  use npop from clusters lfc short for
  ## log2 fold change
  lfc_level <- 0.6
  deinfo <- 1:npop %>% map(function(ipop) {
    symsim_get_degenes(
      true_counts_res = true_counts, popA = ipop,
      popB = setdiff(1:npop, ipop)
    )
  })
  degenes <- deinfo %>% map(function(i) {
    (i$nDiffEVF > 0) & (abs(i$logFC_theoretical) > lfc_level)
  })
  return(list(
    deinfo = deinfo, degenes = degenes, kcc = cc, npop = npop,
    tpm = tpm_gbc, logtpm = logtpm_gbc, tptth = tptth_gbc,
    logtptth = logtptth_gbc, true_counts = true_counts, obs_nonumi = obs_nonumi,
    obs_umi = obs_umi
  ))
}
## gt the theoretic expression values maybe not used here.
get_mean_counts <- function(kinetic_params) {
  kon <- kinetic_params[[1]]
  koff <- kinetic_params[[2]]
  s <- kinetic_params[[3]]
  return(s * kon / (kon + koff))
}


## ** evaluating clusters, DE genes, and imputation mse ***
## clustering aRI kcc short for known cell clusters
get_clustering_ari <- function(logtpm, kcc, k = 2) {
  npop <- nlevels(as.factor(kcc))
  if (is.null(k)) {
    k <- npop - 1
  }
  slogtpm <- t(scale(t(logtpm)))
  slogtpm[is.nan(slogtpm)] <- 0
  s <- svd(slogtpm)
  c <- kmeans(s$v[, 1:k], npop, iter.max = 80, nstart = 300)
  return(mclust::adjustedRandIndex(c$cluster, kcc))
}

## *** AUC of DE genes analysis
get_de_auc <- function(logtpm, kcc, degenes, whichtest = "t") {
  zhu_ttest <- function(x, m, celltype, whichtest) {
    if (sd(x) < 1e-06) {
      1
    } else {
      if (whichtest == "t") {
        t.test(x ~ celltype == m)$p.value
      } else if (whichtest == "wilcox") {
        wilcox.test(x ~ celltype == m)$p.value
      } else {
        stop(paste0("Wrong test name: ", whichtest))
      }
    }
  }
  npop <- length(degenes)
  ## return a vector size of npop
  1:npop %>% map_dbl(function(i) {
    x <- apply(logtpm, 1, zhu_ttest, m = i, celltype = kcc, whichtest = whichtest)
    adjx <- p.adjust(x, method = "fdr")
    caTools::colAUC(X = adjx, y = degenes[[i]])
  })
}

## *** mean square error (mse) of imputation This could be
## done, but seems not that make sense, so skip it now.

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

## need original count data.
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
    imp_count <- read.csv(file = rsf, row.names = 1)
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

## pltfm short for platform
scvi_impute <- function(myseed = 1, pltfm = "umi", ...) {
  if (pltfm == "umi") {
    scale_level <- 10000
  } else if (pltfm == "nonumi") {
    scale_level <- 1e+06
  } else {
    stop("pltfm should be umi or nonumi.")
  }

  fnm <- stringr::str_c("sim", p$ncell, p$ngene, myseed, pltfm, ".csv", sep = "_")
  cnt <- t(read.csv(here::here("simutool", "jobs", "scvi_result", fnm),
                  header=FALSE))
  tpm <- scale(cnt, center = F, scale = colSums(cnt)) * scale_level
  return(log(tpm + 1))
}

## * run method and evaluate it.
main_impute <- function() {
  de_auc_umi <- matrix(0, p$nrep, p$npop)
  de_auc_nonumi <- matrix(0, p$nrep, p$npop)
  npc <- length(p$pcs)
  clust_ari_umi <- matrix(0, p$nrep, npc)
  clust_ari_nonumi <- matrix(0, p$nrep, npc)
  for (i in 1:p$nrep) {
    d <- load_symsim_data(i)
    ## umi experiment
    logtpm_umi <- do.call(
      what = str_c(p$method, "_impute"),
      args = list(
        method = p$method, tpm = d$tptth, logtpm = d$logtptth,
        k = p$k, m = p$m, ncores = p$ncores, obs_count = d$obs_umi$counts,
        scale_level = p$scale_level_umi
      )
    )
    de_auc_umi[i, ] <- get_de_auc(logtpm_umi, d$kcc, d$degenes, whichtest = p$detest)
    for (j in 1:npc) {
      clust_ari_umi[i, j] <- get_clustering_ari(
        logtpm_umi,
        d$kcc, p$pcs[j]
      )
    }
    ## nonumi experiment
    logtpm_nonumi <- do.call(
      what = str_c(p$method, "_impute"),
      args = list(
        method = p$method, tpm = d$tpm, logtpm = d$logtpm,
        k = p$k, m = p$m, ncores = p$ncores, obs_count = d$obs_nonumi$counts,
        scale_level = p$scale_level_nonumi
      )
    )
    de_auc_nonumi[i, ] <- get_de_auc(logtpm_nonumi, d$kcc, d$degenes, whichtest = p$detest)
    for (j in 1:npc) {
      clust_ari_nonumi[i, j] <- get_clustering_ari(
        logtpm_nonumi,
        d$kcc, p$pcs[j]
      )
    }
  }
  result <- list(
    de_auc_umi = de_auc_umi,
    de_auc_nonumi = de_auc_nonumi,
    clust_ari_umi = clust_ari_umi,
    clust_ari_nonumi = clust_ari_nonumi
  )
  saveRDS(object = result, file = paste0(
    "result/",
    str_c(p$method, p$k, p$m, ".rds", sep = "_")
  ))
}

## * script main interface
registerDoParallel(cores = p$ncores)
main_impute()
