# simulation according to scIMpute
library(scImpute)
library(Rmagic)
library(foreach)
library(doParallel)
library(SIMPLE)
library(matrixStats)

source("Simulation/utils.R")

load("../../chu_data/data.rdat")
trueDE <- read.csv("../../chu_data/bulk_diff.csv", row.names = 1)

K0 = 10
M0 = 7
pm = 0.6 #0.7, 0.6; opt2: dropP = 0.2, pm = 0.4; dropP = 0.4, pm = 0.6
dropP  = 0.3 #0.2, 0.4
frac = 0.75
opt = 1

run = 5
cluster_result <- matrix(0, run, 16)
mse <- matrix(0, run, 7)
auc <- matrix(0, run, 6)
mse_cor <- matrix(0, run, 12) # for true_correlation =0  & !=0

colnames(cluster_result) <- paste(rep(c('Y2', 'Y', 'scimpclu_clus','scimpclu_impt', 'scimpclu_bulk_clus','scimpclu_bulk_impt','scimpute', 'magic'), each = 2), c(1,2), sep="-")
colnames(mse) =c('Y2','scimpclu_mean','scimpclu_impt', 'scimpclu_bulk_mean','scimpclu_bulk_impt','scimpute', 'magic')
colnames(auc) = c('Y2', 'Y', 'scimpclu_impt', 'scimpclu_bulk_impt','scimpute', 'magic')
colnames(mse_cor) = paste(rep(c('Y2', 'Y', 'scimpclu', 'scimpclu_bulk','scimpute', 'magic'), each=2), c(1,2), sep="-")

cl = as.numeric(factor(celltype))

for( r in 1:run)
{
  # add more dropouts to original data and subsample cells
  cell_ind = sample(1:ncol(dat_norm), ncol(dat_norm) * frac)
  Y2 = dat_norm[, cell_ind]
  celltype_true = celltype[cell_ind]
  cl_true = cl[cell_ind]
  # add dropout
  G = nrow(Y2)
  n = ncol(Y2)
  #cl = as.numeric(factor(celltype_true))
  if(opt==1)
  {
    Z = matrix(rbinom(G*n,1,exp(-dropP*rowMeans(dat_norm))), nrow=G) #0.3 * exp,
  }else{
    Z = matrix(rbinom(G*n,1,dropP), nrow=G)
  }
  Y2[Z==1] = 0

  # random sample 3000 genes
  ind = sample(which(rowSums(Y2 != 0) > 50), 3000)
  Y2 = Y2[ind, ]
  Y = dat_norm[ind, cell_ind]
  Z2 = Z[ind, ]
  bulk_mean = bulk_norm_mean[ind,]

  # no use?
  prop = xtabs(~celltype_true)/n
  bulk = bulk_mean %*% prop

  # get bulk
  DE = matrix(0, 3000, M0)
  DE[which(abs(trueDE[ind,]) < 1e-6)] = 1

  cluster_result[r,1:2] =getCluster(Y2, celltype_true, Ks = 20, M0 = M0)[[1]]
  cluster_result[r,3:4] =getCluster(Y, celltype_true, Ks = 20, M0 = M0)[[1]]
  mse[r,1] <- mean(Y[Z2==1]^2) #

  Y2_scale = t(scale(t(Y2)))
  s = svd(Y2_scale)
  km0 <- kmeans(t(Y2_scale)%*% s$u[,1:20], M0, iter.max = 80, nstart = 300) #s$v[,1:K]
  #km0 <- kmeans(s$v[,1:20], M0, iter.max = 80, nstart = 300)
  clus = km0$cluster


  #### scimpclu ####
  registerDoParallel(cores = 7)
  result <- scimpclu(Y2, K0, M0, clus = clus, K = 20, iter= 10, est_z = 5, impt_it = 1, max_lambda=T, est_lam = 1, penl = 1, sigma0 = 100, p_min = pm, cutoff=0.01,verbose=T, min_gene = 500, num_mc=5) #as.numeric(factor(celltype_true)) NULL, est_z = 1 clus
  resultY = result$Y

  cluster_result[r,5:6] = mclust::adjustedRandIndex(apply(result$z,1, which.max), celltype_true)
  imptY = result$Yimp0
  imptY[imptY < 0] = 0
  mse[r,2] = mean((imptY[Z2==1] - Y[Z2==1])^2) # should use this

  # plot(rowMeans(Y==0), sapply(1:nrow(Y),function(i) {
  #   ii = which(Z2[i, ]==1)
  #   mean((imptY[i, ii] - Y[i,ii])^2)
  #   }), ylab="mse")

  result1 = do_impute(Y2, result$Y, result$beta, result$lambda, result$sigma, result$mu, result$pi, result$geneM, result$geneSd, result$initclus, mcmc=50, burnin = 5, verbose=F, pg = result$pg, cutoff = 0.01) # only for AUC and clustering

  cluster_result[r,7:8] = getCluster(result1$impt, celltype_true, Ks = 20, M0 = M0)[[1]]# should use this
  mse[r,3] = mean((result1$impt[Z2==1] - Y[Z2==1])^2)

  resultB <- scimpclu_bulk(Y2, K0, bulk, M0, celltype=rep(1, n), clus = clus, K = 20, iter= 10, est_z = 1, impt_it = 1, max_lambda=T, est_lam = 2, penl =1, sigma0 = 100, p_min = pm, min_gene = 300,cutoff=0.01,verbose=T, num_mc=5) # NULL
  resultYB = resultB$Y

  cluster_result[r,7:8] = getCluster(resultB$impt, celltype_true, Ks = 20, M0 = M0)[[1]]# should use this

  imptY = resultB$Yimp0
  imptY[imptY < 0] = 0
  mse[r,2] = mean((imptY[Z2==1] - Y[Z2==1])^2) # should use this

  result2 = do_impute(Y2, resultB$Y, resultB$beta, resultB$lambda, resultB$sigma, resultB$mu, resultB$pi, resultB$geneM, resultB$geneSd, rep(1, n), mcmc=50, burnin = 5, verbose=F, pg = resultB$pg, cutoff = 0.01)

  cluster_result[r,11:12] = getCluster(result2$impt, celltype_true, Ks = 20, M0 = M0)[[1]]
  mse[r,5] = mean((result2$impt[Y2==0] - Y[Y2==0])^2)


  #### method 3 ####
  sci_result <- scimpute("", infile = "csv", outfile = "csv", type = "count",
                         "../scimpute/", labeled = FALSE, drop_thre = pm, Kcluster = 7,
                         labels = NULL, genelen = NULL, ncores =6, count_lnorm = Y2+log10(1.001))

  # plot(rowMeans(Y==0), sapply(1:nrow(Y),function(i) {
  #   ii = which(Z2[i, ]==1)
  #   mean((sci_result[[1]][i, ii] - Y2[i,ii])^2)
  # }), ylab="mse")


  cluster_result[r,13:14] = getCluster(sci_result[[1]], celltype_true, Ks = 20, M0 = M0)[[1]]
  mse[r,6] = mean((sci_result[[1]][Z2==1] - Y[Z2==1])^2)


  # save(Y2, Z2, sci_result, file="chu_data/scimpute2.rdat")
  # save(Y2, sci_result, file="chu_data/scimpute_part_cells_K6.rdat")
  # save(Y2, Y, Z2, cell_ind, result, result2, file="chu_data/fa_impute_part_cells_K6_bulk.rdat")

  #sci_result <- read.csv("chu_data/scimpute_count.csv", row.names = 1)

  #### method 2 ####
  MAGIC <- magic(t(Y2), genes="all_genes", t="auto")
  cluster_result[r,15:16] = getCluster(t(MAGIC$result), celltype_true, Ks = 20, M0 = M0)[[1]]
  mse[r,7] = mean((t(MAGIC$result)[Y2==0] - Y[Y2==0])^2)

  # auc
  auc_temp = foreach(m=1:M0, .combine = rbind) %dopar%
  {
    res0 <- apply(Y2, 1, myttest, m = m, celltype=cl_true)
    #   # if(sd(x) <1e-6) {
    #   #   1
    #   # }else{
    #   #   t.test(x~celltype_true==m)$p.value
    #   # }
    #   # }) #, alternative="l"
    res <- apply(Y, 1, myttest, m = m, celltype=cl_true)
    res2 <- apply(result1$impt, 1, myttest, m = m, celltype=cl_true)
    #res3 <- apply(result2$impt, 1, myttest, m = m, celltype=cl_true)
    res4 <- apply(sci_result[[1]], 1, myttest, m = m, celltype=cl_true)
    #res5 <- apply(MAGIC$result, 2, myttest, m = m, celltype=cl_true)


    c(caTools::colAUC(res0, DE[,m]),
    caTools::colAUC(res, DE[,m]),
    caTools::colAUC(res2, DE[,m]),
    caTools::colAUC(res3, DE[,m]),
    caTools::colAUC(res4, DE[,m])
    ,caTools::colAUC(res5, DE[,m])
    )

  }

  auc[r,] = colSums(auc_temp)


  # gene-gene correlation
  # only consider within cluster correlation, separate no from with correlation
  for(m in 1:M0)
  {
    cors = cor(t(dat_norm[ind,cl==m]))
    cors[is.na(cors)]=0
    mse_cor[r, 1:2] = mse_cor[r, 1:2] + mycor2(Y2, m, cl_true, cors)
    # #hist(org_cors[cors==0], breaks = 50)
    # #hist(org_cors[BB!=0] - cors[BB!=0])
    #
    #mse_cor[r, 3:4] = mse_cor[r, 3:4] + mycor(Y, m, celltype_true, cors)
    mse_cor[r, 5:6] = mse_cor[r, 5:6] + mycor2(resultY, m, cl_true, cors)
    #mse_cor[r, 7:8] = mse_cor[r, 7:8] + mycor2(resultYB, m, cl_true, cors)
    mse_cor[r, 9:10] = mse_cor[r, 9:10] + mycor2(sci_result[[1]], m, cl_true, cors)
    #mse_cor[r, 11:12] = mse_cor[r, 11:12] + mycor2(t(MAGIC$result), m, cl_true, cors)

  }


}








# heatmap of orignial and imputed dataset
# 3 levels of imputation
imputeY = result$Y*result$geneSd + result$geneM # Y, impute for 0
imputeY = t(result$Yimp0)*result$posSd + result$posM # E(Y), impute for all data
imputeY = result2$impt # E(Y) impute for 0

imputeY = t(data_MAGIC$result) #

rownames(imputeY) = 1:nrow(imputeY)
colnames(imputeY) = 1:ncol(imputeY)


sel_gene0 = which(rowSums(impute_DE0) > 0)
sel_gene1 = which(rowSums(impute_DE1[,3:4]) > 0)
sel_gene2 = which(rowSums(impute_DE2) > 0)
sel_gene3 = which(rowSums(impute_DE3) > 0)

pdf("chu_data/DE_genes_venn.pdf")
gplots::venn(list("bulk" = sel_gene0, "faimpute" = sel_gene1,  "org" = sel_gene3))#"scimpute" = sel_gene2,
dev.off()

#figures2 <- heatmap(resultB$Yimp0[sel_gene0,], factor(clus), factor(celltype), grp_ord = NULL, scale=T,cell.ord ="celltype") # use most imputed data to order gene!, ,
figures <- heatmap(result1$impt[sel_gene1, ], factor(clus), factor(celltype), grp_ord = NULL, scale=T,cell.ord ="celltype")#gene_ord = figures2[[2]], ,disp.max = 3

figures0 <- heatmap(Y2[sel_gene1, ], factor(celltype), factor(celltype), grp_ord = NULL, scale=T,cell.ord ="celltype",gene_ord = figures[[2]])#, disp.max = 3

figures2 <- heatmap(sci_result[[1]][sel_gene1,], factor(clus2), factor(celltype), grp_ord = NULL, scale=T,cell.ord ="celltype",gene_ord = figures[[2]]) # use most imputed data to order gene!, ,



plot_grid(figures[[1]][[2]],figures[[1]][[1]], nrow = 2, align = "v",rel_heights=c(1,7)) #figures[[3]],
pdf("scimpute_DE_heatmap_0713.pdf") #chu_data/
plot_grid(figures2[[1]][[2]],figures2[[1]][[1]], nrow = 2, align = "v",rel_heights=c(1,7))
dev.off()
plot(figures0[[1]][[1]])


load("Chu/Chu_example1_dropP-0.4_frac-0.5_viper_repeat-1_option-0.RData")
load("Chu/Chu_example1_dropP-0.4_frac-0.5_viper.RData")

viper_result_imputed_log = log(exp(viper_result$imputed_log) - 0.1 + 1)
