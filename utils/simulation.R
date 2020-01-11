# simulations
library(scImpute)
library(Rmagic)
library(foreach)
library(doParallel)
library(SIMPLEs)
library(msm)
library(ggplot2)
library(reshape2)
library(matrixStats)

# make sure current work directory
# TODO: use here
source("utils/utils.R")

args = commandArgs(trailingOnly = TRUE)

#### step1: impute and estimate B ####
run = 10
cluster_result <- matrix(0, run, 16)
mse <- matrix(0, run, 7)
auc <- matrix(0, run, 6)
mse_cor <- matrix(0, run, 12)  # for true_correlation =0  & !=0


colnames(cluster_result) <- paste(rep(c("Y2", "Y", "SIMPLE_clus", "SIMPLE_impt", 
    "SIMPLE_bulk_clus", "SIMPLE_bulk_impt", "scimpute", "magic"), each = 2), c(1, 
    2), sep = "-")
colnames(mse) = c("Y2", "SIMPLE_mean", "SIMPLE_impt", "SIMPLE_bulk_mean", "SIMPLE_bulk_impt", 
    "scimpute", "magic")
colnames(auc) = c("Y2", "Y", "SIMPLE_impt", "SIMPLE_bulk_impt", "scimpute", "magic")
colnames(mse_cor) = paste(rep(c("Y2", "Y", "SIMPLE", "SIMPLE_bulk", "scimpute", "magic"), 
    each = 2), c(1, 2), sep = "-")


K0 = as.numeric(args[1])
n = 300
M0 = 3
pm = as.numeric(args[2])  #0.5, 0.4
for (r in 1:run) {
    set.seed(r)
    print(paste("run:", r))
    
    #### simualtion 1 #### simu_data = simulation_bulk(n, S0 = 20, K = 0, MC=M0,
    #### block_size = 50, indepG = 1000 , verbose=T)
    
    
    #### simualtion 2 ####
    simu_data = simulation_bulk(n, S0 = 20, K = 6, MC = M0, block_size = 32, indepG = 1000 - 
        32 * 6, verbose = F, overlap = 0)
    # simu_data = simulation_bulk(n, S0 = 20, K = 3, MC=M0, block_size = 100, indepG
    # = 1000 - 190, verbose=F, overlap=55)
    
    
    Y2 = simu_data$Y2
    Y = simu_data$Y
    label_true = simu_data$S_label
    celltype_true = simu_data$Z
    cluster_result[r, 1:2] = getCluster(Y2, celltype_true, Ks = 20, M0 = M0)[[1]]
    cluster_result[r, 3:4] = getCluster(Y, celltype_true, Ks = 20, M0 = M0)[[1]]
    mse[r, 1] <- mean(Y[Y2 == 0]^2)
    
    #### SIMPLE ####
    registerDoParallel(cores = 6)
    result <- SIMPLE(Y2, K0, M0, clus = NULL, K = 20, iter = 10, est_z = 1, impt_it = 1, 
        max_lambda = T, est_lam = 1, penl = 1, sigma0 = 100, p_min = pm, cutoff = 0.01, 
        verbose = T, min_gene = 200)  #0.4
    resultY = result$Y
    
    
    cluster_result[r, 5:6] = mclust::adjustedRandIndex(apply(result$z, 1, which.max), 
        celltype_true)
    mse[r, 2] = mean((result$Yimp0[Y2 == 0] - Y[Y2 == 0])^2)  # should use this
    
    result1 = do_impute(Y2, result$Y, result$beta, result$lambda, result$sigma, result$mu, 
        result$pi, result$geneM, result$geneSd, result$initclus, mcmc = 50, burnin = 5, 
        verbose = F, pg = result$pg, cutoff = 0.01)  # only for AUC
    
    cluster_result[r, 7:8] = getCluster(result1$impt, celltype_true, Ks = 20, M0 = M0)[[1]]  # should use this
    mse[r, 3] = mean((result1$impt[Y2 == 0] - Y[Y2 == 0])^2)
    
    #### using bulk ####
    result <- SIMPLE_B(Y2, K0, data.frame(simu_data$bulk), M0, celltype = rep(1, 
        n), clus = NULL, K = 20, iter = 10, est_z = 1, impt_it = 1, max_lambda = T, 
        est_lam = 2, penl = 1, sigma0 = 100, p_min = pm, min_gene = 200, cutoff = 0.01, 
        verbose = T)
    resultYB = result$Y
    
    cluster_result[r, 9:10] = mclust::adjustedRandIndex(apply(result$z, 1, which.max), 
        celltype_true)  #getCluster(result$Yimp0, celltype_true, Ks = 20, M0 = M0)[[1]]
    mse[r, 4] = mean((result$Yimp0[Y2 == 0] - Y[Y2 == 0])^2)
    
    result2 = do_impute(Y2, result$Y, result$beta, result$lambda, result$sigma, result$mu, 
        result$pi, result$geneM, result$geneSd, rep(1, n), mcmc = 50, burnin = 5, 
        verbose = F, pg = result$pg, cutoff = 0.01)
    
    cluster_result[r, 11:12] = getCluster(result2$impt, celltype_true, Ks = 20, M0 = M0)[[1]]
    mse[r, 5] = mean((result2$impt[Y2 == 0] - Y[Y2 == 0])^2)
    
    
    #### method 3 ####
    sci_result <- scimpute("", infile = "csv", outfile = "csv", type = "count", "../scimpute/", 
        labeled = FALSE, drop_thre = pm, Kcluster = M0, labels = NULL, genelen = NULL, 
        ncores = 6, count_lnorm = data.frame(Y2 + log10(1.001)))
    
    
    cluster_result[r, 13:14] = getCluster(sci_result[[1]], celltype_true, Ks = 20, 
        M0 = M0)[[1]]
    mse[r, 6] = mean((sci_result[[1]][Y2 == 0] - Y[Y2 == 0])^2)
    
    #### method 4 ####
    MAGIC <- magic(t(Y2), genes = "all_genes", t = "auto")
    cluster_result[r, 15:16] = getCluster(t(MAGIC$result), celltype_true, Ks = 20, 
        M0 = M0)[[1]]
    mse[r, 7] = mean((t(MAGIC$result)[Y2 == 0] - Y[Y2 == 0])^2)
    
    
    # auc
    for (m in 1:M0) {
        res0 <- apply(Y2, 1, myttest, m = m, celltype = celltype_true)
        res <- apply(Y, 1, myttest, m = m, celltype = celltype_true)
        res2 <- apply(result1$impt, 1, myttest, m = m, celltype = celltype_true)
        res3 <- apply(result2$impt, 1, myttest, m = m, celltype = celltype_true)
        res4 <- apply(sci_result[[1]], 1, myttest, m = m, celltype = celltype_true)
        res5 <- apply(MAGIC$result, 2, myttest, m = m, celltype = celltype_true)
        
        
        auc[r, 1] = auc[r, 1] + caTools::colAUC(res0, label_true[, m])
        auc[r, 2] = auc[r, 2] + caTools::colAUC(res, label_true[, m])
        auc[r, 3] = auc[r, 3] + caTools::colAUC(res2, label_true[, m])
        auc[r, 4] = auc[r, 4] + caTools::colAUC(res3, label_true[, m])
        auc[r, 5] = auc[r, 5] + caTools::colAUC(res4, label_true[, m])
        auc[r, 6] = auc[r, 6] + caTools::colAUC(res5, label_true[, m])
        
    }
    
    # only consider within cluster correlation, separate no from with correlation
    BB = simu_data$B %*% t(simu_data$B)
    cors = (BB + diag(simu_data$Sigma^2))/sqrt(diag(BB) + simu_data$Sigma^2)
    cors = t(cors)/sqrt(diag(BB) + simu_data$Sigma^2)
    
    for (m in 1:M0) {
        mse_cor[r, 1:2] = mse_cor[r, 1:2] + mycor(Y2, m, celltype_true, cors)
        mse_cor[r, 3:4] = mse_cor[r, 3:4] + mycor(Y, m, celltype_true, cors)
        mse_cor[r, 5:6] = mse_cor[r, 5:6] + mycor(resultY, m, celltype_true, cors)
        mse_cor[r, 7:8] = mse_cor[r, 7:8] + mycor(resultYB, m, celltype_true, cors)
        mse_cor[r, 9:10] = mse_cor[r, 9:10] + mycor(sci_result[[1]], m, celltype_true, 
            cors)
        mse_cor[r, 11:12] = mse_cor[r, 11:12] + mycor(t(MAGIC$result), m, celltype_true, 
            cors)
        
    }
    
    
}

save(cluster_result, mse, auc, mse_cor, file = paste("SIMPLE/Simulation/simulation_K6_K0", 
    K0, "pm", pm, "lowdrop.rdat", sep = "_"))

# summarizing the results
round(colMeans(cluster_result[, seq(1, 16, by = 2)]), 2)
round(colSds(cluster_result[, seq(1, 16, by = 2)])/sqrt(10), 3)

round(colMeans(mse), 2)
round(colSds(mse)/sqrt(10), 2)

round(colMeans(auc/M0), 2)
round(colSds(auc/M0)/sqrt(10), 2)

round(colMeans(mse_cor[, seq(1, 12, by = 2)]/M0), 3)
round(colSds(mse_cor[, seq(1, 12, by = 2)]/M0)/sqrt(10), 3)

round(colMeans(mse_cor[, seq(2, 12, by = 2)]/M0), 3)
round(colSds(mse_cor[, seq(2, 12, by = 2)]/M0)/sqrt(10), 3)




