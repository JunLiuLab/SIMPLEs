## ref simutool/eval_simulation_singlejob.R
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
  make_option(c("-t", "--test"),
    action = "store",
    type = "character", default = "t"
  ),
  make_option(c("-n", "--ncores"),
    action = "store",
    type = "integer", default = 1
  ),
  make_option(c("--thres"),
    action = "store",
    type = "double", default = 2.0
  )
)

args <- option_list %>%
  OptionParser(option_list = .) %>%
  parse_args()


## * load data
sccnt <- readRDS(here::here("10xGenomics", "scRNAseq", "sc_10x_5cl_rgmm_cnt.rds"))
bulkinfo <- readRDS(here::here("10xGenomics", "bulk_RNAseq", "DESeq2_result.RDS"))
cell_cluster <- gsub(".*:", "", colnames(sccnt))
sc_genes <- rownames(sccnt)
bk_genes <- bulkinfo$gmap
## ** join and align:  scgnes, bk_genes
common_genes <- intersect(sc_genes, bk_genes)

## rr short for rotate rows
sccnt_rr <- sccnt[common_genes, ]

## ** save scrnaseq count for scvi
for_scvi_cnt <- function() {
  write.csv(sccnt_rr, file = here::here(
    "10xgenomics", "scrnaseq",
    "sc_10x_5cl_forimput_cnt.csv"
  ))
}
## ** get bulk differential genes
get_bulk_degs <- function(deseq_res,
                          reduce_row_names = common_genes,
                          log2fc_thres = args$thres) {
  ## reduce rows by gene symbols.
  t <- deseq_res[reduce_row_names, ]
  r <- ((abs(t$log2fc) > log2fc_thres) & (t$padj < 0.05))
  print("differential gene number is ")
  print(length(which(r == TRUE)))
  return(r)
}

h1975va549 <- get_bulk_degs(bulkinfo$h1975va549)
h2228va549 <- get_bulk_degs(bulkinfo$h2228va549)
h838va549 <- get_bulk_degs(bulkinfo$h838va549)
hcc827va549 <- get_bulk_degs(bulkinfo$hcc827va549)


## * set global parameters
scale_level <- 10000
pcs <- seq(2, 50, 2)

## * utils
## ** seq data transfromation
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

## ** measure clustering and DE analysis
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

# group: logic vector, length of row nm of logtpm
get_de_auc <- function(logtpm, mytest = "t") {
  zhu_ttest <- function(x, group) {
    if (sd(x) < 1e-06) {
      1
    } else {
      if (mytest == "t") {
        t.test(x ~ group)$p.value
      } else if (mytest == "wilcox") {
        wilcox.test(x ~ group)$p.value
      } else {
        stop(paste0("Wrong test name: ", mytest))
      }
    }
  }

  two_group <- function(testg = "H1975", degenes, control = "A549") {
    kept_cols <- cell_cluster %in% c(testg, control)
    logtpm_2grp <- logtpm[, kept_cols]
    group_2grp <- gsub(".*:", "", colnames(logtpm_2grp)) %in% testg
    x <- apply(logtpm_2grp, 1, zhu_ttest, group = group_2grp)
    adjx <- p.adjust(x, method = "fdr")
    caTools::colAUC(X = adjx, y = degenes)
  }

  de_auc_umi <- rep(0, 4)
  de_auc_umi[1] <- two_group("H1975", h1975va549)
  de_auc_umi[2] <- two_group("H2228", h2228va549)
  de_auc_umi[3] <- two_group("H838", h838va549)
  de_auc_umi[4] <- two_group("HCC827", hcc827va549)
  return(de_auc_umi)
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

## pltfm short for platform
## scvi need pytorch and python env.
## get the results from a seperated job.
scvi_impute <- function(...) {
  cnt <- t(read.csv(here::here("10xGenomics", "impt", "scvi_impt.csv"),
    header = FALSE
  ))
  tpm <- scale(cnt, center = F, scale = colSums(cnt)) * scale_level
  a <- log(tpm + 1)
  rownames(a) <- rownames(sccnt_rr)
  colnames(a) <- colnames(sccnt_rr)
  return(a)
}

## * run method and evaluate it.
main_impute <- function() {
  ## ** cnt transformation
  tpm_umi <- get_tpm(sccnt_rr, scale_level)
  logtpm_umi <- log_tpm_rseq(tpm_umi)
  ## ** parameters
  method <- args$method
  k <- args$K0
  m <- args$M0
  t <- args$test
  ## ** imputation
  logtpm_umi <- do.call(
    what = str_c(args$method, "_impute"),
    args = list(
      method = method, tpm = tpm_umi, logtpm = logtpm_umi,
      k = k, m = m, ncores = args$ncores, obs_count = sccnt_rr,
      scale_level = scale_level
    )
  )

  ## ** save impute results
  saveRDS(logtpm_umi, file = here::here(
    "10xGenomics", "impt",
    str_c(method, k, m, "impt.rds", sep = "_")
  ))

  de_auc_umi <- get_de_auc(logtpm_umi, mytest = t)
  print(de_auc_umi)

  cl_ari_umi <- pcs %>% map_dbl(function(pc) {
    get_clustering_ari(logtpm_umi, cell_cluster, pc)
  })
  print(cl_ari_umi)

  result <- list(
    de_auc_umi = de_auc_umi,
    cl_ari_umi = cl_ari_umi
  )
  saveRDS(object = result, file = here::here(
    "10xGenomics", "result",
    str_c(method, k, m, t, ".rds", sep = "_")
  ))
}

de_after_impute <- function() {
  ## ** parameters
  method <- args$method
  k <- args$K0
  m <- args$M0
  t <- args$test
  c <- args$thres

  ## ** load impute results
  logtpm_umi <- readRDS(here::here(
    "10xGenomics", "impt",
    str_c(method, k, m, "impt.rds", sep = "_")
    ))
  ## double check the names.
  ## some results might have no names and only data.
  rownames(logtpm_umi) <- rownames(sccnt_rr)
  colnames(logtpm_umi) <- colnames(sccnt_rr)

  de_auc_umi <- get_de_auc(logtpm_umi, mytest = t)
  print(de_auc_umi)

  saveRDS(object = de_auc_umi, file = here::here(
    "10xGenomics", "result",
    str_c(method, k, m, t, c, ".rds", sep = "_")
  ))
}

## * script main interface
registerDoParallel(cores = args$ncores)
# after impute, main_impute()
de_after_impute()
