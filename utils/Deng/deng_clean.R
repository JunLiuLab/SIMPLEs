library(scImpute)
library(Rmagic)
library(foreach)
library(doParallel)
library(SIMPLE)
library(matrixStats)

source("Simulation/utils.R")

load("Datasets/deng_dat.rdat")

args = commandArgs(trailingOnly=TRUE)

#### imputation ####
K0=as.numeric(args[1])
M0=6
pm= as.numeric(args[2])
penl = as.numeric(args[3])

# Y2_scale = t(scale(t(Y2)))
# s = svd(Y2_scale)
# km0 <- kmeans(t(Y2_scale)%*% s$u[,1:5], 6, iter.max = 80, nstart = 300)
# #km0 <- kmeans(s$v[,1:10], M0+1, iter.max = 80, nstart = 300)
# clus = km0$cluster
# mclust::adjustedRandIndex(clus, celltype_true)
# xtabs(~clus + celltype_true)
# celltype_true = celltype_true[clus!=4]

registerDoParallel(cores = 8)
result <- scimpclu(Y2, K0, M0, clus = clus, K = 10, iter= 20, est_z = 10, impt_it = 1, max_lambda=T, est_lam = 1, penl = penl, sigma0 = 100, p_min = pm, cutoff=0.01,verbose=F, min_gene = 3000, num_mc=3)

xtabs(~apply(result$z, 1,which.max) + celltype_true)

result1 = do_impute(Y2, result$Y, result$beta, result$lambda, result$sigma, result$mu, result$pi, result$geneM, result$geneSd, result$initclus, mcmc=90, burnin = 10, verbose=T, pg = result$pg, cutoff = 0.01)

save(result, result1, uncer, file="SIMPLE/Deng/deng_SIMPLE_0712-8-0.8-0.2-multi.rdat")

# get uncertainty score for each cell
uncer = result1$consensus_cluster
uncer[uncer > 0.5] = 1 - uncer[uncer > 0.5]
uncer = rowMeans(uncer)

getCluster(result1$impt, celltype_true, Ks = c(5, 10,15), M0 = M0)
getCluster(result1$impt, cl2, Ks = c(10:15), M0 = 8)

save(result, result1, K0, M0, pm, file=paste("deng_SIMPLE_0712",K0, pm,penl,"initclus.rdat", sep = "-"))

clus2 = apply(result$z, 1, which.max)
tsneplot(t(result1$impt), factor(cl2), uncer, scale=T, center=T, K=11, file="deng_simple_tsne-cl2-multi.pdf") #


