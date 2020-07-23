## ---- echo = FALSE------------------------------------------------------------
  knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## -----------------------------------------------------------------------------
library(SIMPLEs)
library(doParallel)
library(foreach)
# SIMPLEs may need parallel to speed up.
registerDoParallel(cores = 8)

## ----eval=FALSE---------------------------------------------------------------
#  library(SIMPLEs)
#  # set cluster number, default is 1
#  M0 <- 1
#  # set latent module dimenstion, default is 10.
#  K0 <- 10
#  result <- SIMPLE(dat = data, K0=K0, M0=M0)
#  # result$impt is the imputated result for the data.

## -----------------------------------------------------------------------------
simulation_bulk <- function(n = 300, S0 = 10, K = 3, MC = 2, block_size = 50, overlap = 15, indepG = 30, dropout = 0.3) {
    B = NULL
    W = NULL
    Lambda = NULL
    # sample Y from Factor model, Y = BW + E, W: K*n, Y: G*n, E: G*n
    Z = rmultinom(n, 1, rep(1/MC, MC))
    Z = apply(Z, 2, which.max)
    Z = sort(Z)
    
    if (K > 0) {
        # generate data
        G = block_size + (K - 1) * (block_size - overlap) + indepG  #120
    } else {
        G = indepG
    }
    Sigma = rgamma(G, 2, 1/0.3)  #rnorm(G, 0.6, 0.1)
    Mu = exp(rnorm(G, mean = 0.5, sd = 0.5)) %*% t(rep(1, MC))
    
    act_ind = sample(1:G, S0 * MC)
    label0 = matrix(0, G, MC)
    # specify mean for each cluster
    for (m in 1:MC) {
        ss = ((m - 1) * S0 + 1):(m * S0)
        Mu[act_ind[ss], m] = Mu[act_ind[ss], m] * sample(c(0.2, 0.5, 1.2, 1.5, 2), S0, replace = T)
        label0[act_ind[ss], m] = 1
    }
    if (K > 0) {
        # Factor Loading
        Gamma <- matrix(0, G, K)
        Gamma[1:block_size, 1] <- 1
        if (K > 1) {
            for (i in 1:K) {
                Gamma[((block_size - overlap) * (i - 1) + 1):
                        (block_size + (block_size - overlap) * (i - 1)), i] <- 1
            }
        }
        B = Gamma/4  #sqrt(block_size) # eigenvalue of BB^T = block_size * B^2
        # specify variance for each cluster
        Lambda = list()
        for (m in 1:MC) {
            Lambda[[m]] = diag(1, K, K)
        }
        # matrix(rnorm(K*n), K, n) # K * n
        W = sapply(Z, function(z) mixtools::rmvnorm(1, rep(0, K), Lambda[[z]]))
        E = matrix(rnorm(G * n), nrow = G) * Sigma
        Y = Mu[, Z] + B %*% W + E
    } else {
        Y = matrix(rnorm(G * n, Mu[, Z], Sigma), nrow = G)
    }
    # add dropout
    Y2 = Y
    dZ = matrix(rbinom(G * n, 1, exp(-dropout * rowMeans(Mu)^2)), nrow = G)  # 0.3
    Y2[dZ == 1] = 0
    Y2[Y2 < 0] = 0
    ind = which(rowSums(Y2 != 0) > 4)
    Y2 = Y2[ind, ]
    Y = Y[ind, ]
    label0 = label0[ind, ]
    B = B[ind, ]
    Mu = Mu[ind, ]
    Sigma = Sigma[ind]
    return(list(Y2 = Y2, Y = Y, B = B, W = W, Mu = Mu, 
                Lambda = Lambda, Sigma = Sigma, Z = Z, bulk = rowMeans(Mu[, Z]), S_label = label0))
}

## ---- eval=FALSE--------------------------------------------------------------
#  n_cell <- 300
#  simu_data <- simulation_bulk(n=n_cell, S0=20, K=3, MC=3, block_size=100, indepG=1000-190,
#                               overlap=55, dropout=0.3)
#  celltype_true = simu_data$Z

## ---- eval=FALSE--------------------------------------------------------------
#  simple_res <- SIMPLE(simu_data$Y2, K0=3, M0=3,p_min=0.6,max_lambda = T,cutoff=0.01)
#  # if we have the bulk data
#  simpleb_res <- SIMPLE_B(simu_data$Y2, K0=3, bulk=data.frame(simu_data$bulk),
#                          celltype=rep(1,n_cell), M0=3,p_min=0.6,max_lambda = T,cutoff=0.01)

## ---- eval=FALSE--------------------------------------------------------------
#  # get the cluster infered from SIMPLEs
#  simple_infered_cluster <- apply(simple_res$z, 1, which.max)
#  # if we have use the SIMPLE_B
#  simpleb_infered_cluster <- apply(simpleb_res$z, 1, which.max)

## ---- eval=FALSE--------------------------------------------------------------
#  # re-do the clustering. Only use the result from SIMPlE as an example.
#  scaled_data <- t(scale(t(simple_res$impt)))
#  s <- svd(scaled_data)
#  km <- kmeans(t(scaled_data) %*% s$u[, 1:20], 3, iter.max=80, nstart=300)

## ---- eval=FALSE--------------------------------------------------------------
#  library(mclust)
#  # between 0 and 1; the larger, the better
#  cluster_score <- mclust::adjustedRandIndex(km$cluster, celltype_true)
#  cluster_score1 <- mclust::adjustedRandIndex(simple_infered_cluster, celltype_true)

## ---- eval=FALSE--------------------------------------------------------------
#  library(Rtsne)
#  library(ggplot2)
#  # for reproducible, we set the seed to tsne.
#  set.seed(0)
#  tsne_res <- Rtsne::Rtsne(t(simple_res$impt), pca_scale=T, pca_center=T, initial_dims=20,
#                           pca = T)
#  partial_tsne <- data.frame(cbind(tsne_res$Y[, 1:2]), celltype_true)
#  colnames(partial_tsne) <- c("TSNE1", "TSNE2", "Type")
#  p <- ggplot(partial_tsne, aes(x = TSNE1, y=TSNE2, color=Type)) + geom_point() + theme_bw()
#  plot(p)

## ---- eval=FALSE--------------------------------------------------------------
#  # load chu dataset
#  data(chu)

## ---- eval = FALSE------------------------------------------------------------
#  ncell = ncol(chu$chu_normalized_data)
#  chu_simpleb_res <- SIMPLE_B(chu$chu_normalized_data,
#                              bulk=data.frame(chu$chu_bulk_mean),
#                              celltype=rep(1,ncell), K0=10, M0=7,max_lambda = T)

## ---- eval=FALSE--------------------------------------------------------------
#  scaled_data <- t(scale(t(chu_simpleb_res$impt)))
#  s <- svd(scaled_data)
#  km <- kmeans(t(scaled_data) %*% s$u[, 1:20], 7, iter.max=80, nstart=300)
#  cluster_score = mclust::adjustedRandIndex(km$cluster, chu$chu_cell_type)

## ---- eval=FALSE--------------------------------------------------------------
#  dropP <- 0.3
#  G <- nrow(chu$chu_normalized_data)
#  n <- ncol(chu$chu_normalized_data)
#  Z <- matrix(rbinom(G*n, 1, exp(-dropP * rowMeans(chu$chu_normalized_data))),
#              nrow=G)
#  more_dropout_data <- chu$chu_normalized_data
#  more_dropout_data[Z==1] <- 0

## ---- eval = FALSE------------------------------------------------------------
#  chu_more_dropout_impute <- SIMPLE_B(more_dropout_data, celltype=rep(1,ncell),
#                                      K0=10, M0=7,bulk=data.frame(chu$chu_bulk_mean))
#  # the smaller, the better
#  mse <- mean((chu_more_dropout_impute$impt[Z==1] - chu$chu_normalized_data[Z==1])^2)
#  
#  simpleb_infered_cluster <- apply(chu_more_dropout_impute$z, 1, which.max)
#  mclust::adjustedRandIndex(simpleb_infered_cluster, chu$chu_cell_type)

