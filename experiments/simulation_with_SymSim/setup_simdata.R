## * Simultion tools
library(SymSim)
library(tidyverse)

## * configs
nbatch <- 1
add_batch_effect <- T
batch_effect_size <- 1
ncell <- 300
ngene <- 1000

## ** cell evf settings cell cluster
npopulation <- 5
phyla <- Phyla5()
min_popsize <- 50
i_minpop <- 1

evf_center <- 1 # should always fix as 1
evf_type <- "discrete"
nevf <- 10

# increasing n_de_evfs can also make the clusters more
# distinct as they represent the number of biological
# conditions that cause differences between cell types.
n_de_evf <- 5
vary <- "s" # 'all' or 's'

sigma <- 0.2 # use for seperating cell clusters

## ** gene module setup in two situations
nks <- 3
nkgs <- 100
nkb <- 6
nkgb <- 30

## ** get gene length
data(gene_len_pool)
gene_len <- sample(gene_len_pool, ngene, replace = FALSE)

## ** UMI settings
depth_mean_umi <- 45000
depth_sd_umi <- 4500

alpha_mean_umi <- 0.1

## ** non UMI settings
depth_mean_nonumi <- 1e+05
depth_sd_nonumi <- 10000

alpha_mean_nonumi <- 0.4
## ** utils
my_sim_true <- function(seed = 0) {
  SimulateTrueCounts(
    ncells_total = ncell, min_popsize = min_popsize,
    i_minpop = i_minpop, ngenes = ngene, evf_center = evf_center,
    evf_type = evf_type, nevf = nevf, n_de_evf = n_de_evf,
    impulse = F, vary = vary, Sigma = sigma, phyla = phyla,
    geffect_mean = 0, gene_effects_sd = 1, gene_effect_prob = 0.3,
    bimod = 0, param_realdata = "zeisel.imputed", scale_s = 1,
    prop_hge = 0.015, mean_hge = 5, randseed = seed
  )
}

my_sim_obs <- function(protocol, true_data, add_batch_effect = T,
                       nbatch = 1, batch_effect_size = 1) {
  depth_mean <- ifelse(protocol == "UMI", depth_mean_umi, depth_mean_nonumi)
  depth_sd <- ifelse(protocol == "UMI", depth_sd_umi, depth_sd_nonumi)
  alpha_mean <- ifelse(protocol == "UMI", alpha_mean_umi, alpha_mean_nonumi)
  tmp <- True2ObservedCounts(
    true_counts = true_data[[1]],
    meta_cell = true_data[3], protocol = protocol, alpha_mean = alpha_mean,
    alpha_sd = 0.02, gene_len = gene_len, depth_mean = depth_mean,
    depth_sd = depth_sd, nPCR1 = 14
  )
  if (add_batch_effect) {
    tmp <- DivideBatches(
      observed_counts_res = tmp, nbatch = nbatch,
      batch_effect_size = batch_effect_size
    )
  }
  intc <- tmp
  intc[[1]] <- apply(tmp[[1]], c(1, 2), function(x) {
    ifelse(x > 0, as.integer(x + 1), 0L)
  })
  return(intc)
}

symsim_ptsne <- function(protocol, obs_data, label = "cell_meta.pop") {
  PlotTsne(
    meta = obs_data[[2]], data = log2(obs_data[[1]] +
      1), evf_type = "discrete", n_pc = 20, label = label,
    saving = F, plotname = protocol
  )
}

mysymsim <- function(seed) {
  ## * scRNAseq data simulation ** get true expressions
  true_counts <- my_sim_true(seed)

  ## ** get observed expressions
  obs_umi <- my_sim_obs("UMI", true_counts,
    add_batch_effect = add_batch_effect,
    nbatch = nbatch, batch_effect_size = batch_effect_size
  )
  obs_nonumi <- my_sim_obs("nonUMI", true_counts,
    add_batch_effect = add_batch_effect,
    nbatch = nbatch, batch_effect_size = batch_effect_size
  )
  save(true_counts, obs_umi, obs_nonumi, file = stringr::str_c("symsim_data/sim",
    ncell, ngene, seed, ".RData",
    sep = "_"
  ))
}

main <- function(myrep=10) {
  for (i in 1:myrep) {
    mysymsim(i)
  }
  ## *** obs data view tsne_true <- symsim_ptsne('TrueCounts',
  ## list(true_counts[[1]], true_counts[[3]]), label = 'pop')
  ## tsne_obs_umi <- symsim_ptsne('UMI', obs_umi, label =
  ## 'cell_meta.pop') tsne_obs_nonumi <- symsim_ptsne('nonUMI',
  ## obs_nonumi, label = 'cell_meta.pop')
}

main(myrep=20)
