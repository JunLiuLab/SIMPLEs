library(foreach)
library(doParallel)
# simulate number of clusters
M0 <- 3
# number of cells
n <- 300
# simulation_bulk and getCluster is defined in the util.R under the util directory of the corresponding github repository.
source("utils/utils.R")
simu_data <- simulation_bulk(n = 900, S0 = 20, K = 6, MC = M0, block_size = 32, indepG = 1000 - 32 * 6, verbose = F, overlap = 0)
Y2 <- simu_data$Y2

# number of factors
K <- c(3, 6, 10)
M <- c(1, 3)
# parallel
registerDoParallel(cores = 6)
# estimate the parameters and sample imputed values
results <- selectKM(Y2, K0=K, M0=M, clus = NULL, K = 20, p_min = NULL, max_lambda = T, min_gene = 200, cutoff = 0.01, rel=0.001)
print(sprintf("best M and K: %d, %d", results$mM, results$mK))
result = results$result
# evaluate cluster performance
celltype_true <- simu_data$Z
mclust::adjustedRandIndex(apply(result$z, 1, which.max), celltype_true)
# or redo clustering based on imputed values (sometimes work better for real data)
getCluster(result$impt, celltype_true, Ks = 20, M0 = M0)[[1]]
