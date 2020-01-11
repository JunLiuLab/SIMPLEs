#' Monte Carlo EM algorithm for imputation and clustering
#'
#'Monte Carlo EM algorithm to sample the imputed values, cluster the
#' cells and learn the correlation structure of genes in each cluster.
#'
#' @details
#' Suppose there are G genes and n cells. For each cell cluster, the gene
#'   expression follows \eqn{Y|Z=m~MVN(\mu_m, B\Lambda_m B^T + \Sigma_m)} where
#'   B is a G by K0 matrix, \eqn{\Sigma_m} is a G by G diagonal matrix whose
#'   diagonal entries are specified by \emph{sigma}, and \eqn{\Lambda_m} is a K0
#'   by K0 diagonal matrix whose diagonal entries are specified by
#'   \emph{lambda}. \eqn{P(Z_m) = \pi_m} where \eqn{\pi~Dir(\alpha)}. We remove
#'   the overall mean of each gene before running the algorithm and all the
#'   parameters are estimated based on the normalized gene expression matrix.
#'   The overall mean is returned as \emph{geneM}.
#'
#' @param Y An initial imputed gene expression matrix.
#' @param Y0 Original scRNASeq data matrix.
#' @param pg A matrix for dropout rate of each cell type. Each row is a gene,
#'   each column is the dropout rate of a cell type. The columns should be
#'   ordered as the cell type label in \emph{clus}.
#' @param M0 Number of clusters.
#' @param K0 Number of latent gene modules.
#' @param cutoff The value below cutoff is treated as no expression.
#' @param iter Number of EM steps.
#' @param beta A G by K0 matrix. Initial values for factor loadings (B). See
#'   details.
#' @param sigma A G by M0 matrix. Initial values for the variance of
#'   idiosyncratic noises. Each column is for a cell cluster. See details.
#' @param lambda A M0 by K0 matrix. Initial values for the variances of factors.
#'   Each column is for a cell cluster. See details.
#' @param pi A vector for initial probabilites of cells belong to each cluster.
#' @param z A n by M0 matrix for the probability of each cell belonging to each
#'   cluster. Can be initialized as the one-hot encoding of cluster membership
#'   of cells. If null, z will be updated in the first iteration.
#' @param mu A G by M0 matrix. Initial values for the gene expression mean of
#'   each cluster. Each column is for a cell cluster. If NULL, it will take the
#'   sample mean of cells weighted by the probability in each cluster. See
#'   details.
#' @param celltype A numeric vector for labels of cells in the scRNASeq. Each
#'   cell type has different dropout rate. If input bulk RNASeq data, each cell
#'   type has corresponding mean expression in the bulk RNASeq data. The labels
#'   must start from 1 to the number of types. If NULL, all cells are treated as
#'   a single cell type.
#' @param penl L1 penalty for the factor loadings.
#' @param est_z The iteration starts to update z.
#' @param max_lambda Whether to maximize over lambda.
#' @param est_lam The iteration starts to estimate lambda.
#' @param impt_it The iteration starts to sample new imputed values.
#' @param sigma0 The variance of the prior distribution of \eqn{\mu}.
#' @param pi_alpha The hyperparameter of the prior distribution of \eqn{\pi}.
#'   See details.
#' @param verbose Whether to show some intermediate results. Default = False.
#'
#' @return \code{EM_impute} returns a list of results in the following order.
#' \enumerate{
#'\item{loglik}{The log-likelihood of the imputed gene expression at each iteration.}
#' \item{pi}{Probabilites of cells belong to each cluster.}
#' \item{mu}{Mean expression for each cluster.}
#' \item{sigma}{Variances of idiosyncratic noises for each cluster.}
#' \item{beta}{Factor loadings.}
#' \item{lambda}{Variances of factors for each cluster.}
#' \item{z}{The probability of each cell belonging to each cluster.}
#' \item{Ef}{Conditonal expection the factors for each cluster \eqn{E(f_i|z_i = m)}. A list with length M0, each element in the list is a n by K0 matrix.}
#' \item{Varf}{Conditonal covariance of factors for each cluster \eqn{Var(f_i|z_i = m)}. A list with length M0, each element in the list is a K0 by K0 matrix.}
#' \item{Y}{Last sample of imputed matrix.}
#' \item{geneM}{Overall mean of each gene expression. See details. } %Need to scale Y back to get the imputed expression matrix.
#' \item{geneSd}{Equal to 1 for each gene.}
#' }
#'
#' @author Zhirui Hu, \email{zhiruihu@g.harvard.edu}
#' @author Songpeng Zu, \email{songpengzu@g.harvard.edu}

# if est_z = 1, not use initial values of z
EM_impute <- function(Y, Y0, pg, M0, K0, cutoff, iter, beta, sigma, lambda, pi, z, 
    mu = NULL, celltype = NULL, penl = 1, est_z = 2, max_lambda = T, est_lam = 2, 
    impt_it = 5, sigma0 = 100, pi_alpha = 1, verbose = F, num_mc = 3, lower = -Inf, 
    upper = Inf) {
    
    PI = 3.14159265359
    
    Y[Y > upper] <- upper
    Y[Y < lower] <- lower
    
    n <- ncol(Y)
    G <- nrow(Y)
    
    gene_mean <- rowMeans(Y)
    gene_sd <- rep(1, G)  # rowSds(Y)
    Y <- Y - rowMeans(Y)  # t(scale(t(Y)))
    
    
    if (is.null(mu)) {
        if (M0 > 1 & !is.null(z)) {
            mu <- Y %*% z/n
        } else {
            mu <- matrix(0, G, M0)
        }
    }
    
    if (is.null(celltype)) 
        celltype <- rep(1, n)  # cell type label for bulk and pg
    
    loglik <- matrix(0, n, M0)
    loglik0 <- matrix(0, n, M0)
    ## rownames(loglik0) = colnames(dat)
    tot <- rep(0, iter)
    
    M <- list()  # (BSigma^-1B + Lambda^-1)^-1
    Wm <- list()  # mean of each factor, K0*n
    
    equal <- function(x, Em, nz) {
        sum(x)
    }
    
    fn <- function(x, Em, nz) {
        sum(Em/x + log(x) * (nz - 2))  # -2(alpha-1)
    }
    
    for (it in 1:iter) {
        # E step get expectation of Z
        for (m in 1:M0) {
            beta2 <- beta/sigma[, m]
            M[[m]] <- solve(t(beta) %*% beta2 + diag(1/lambda[m, ]))
            dt <- sum(log(sigma[, m])) + sum(log(lambda[m, ])) - log(det(M[[m]]))  # !!! logdet !=det(log=T)
            loglik0[, m] <- apply(Y, 2, function(x) {
                x2 <- x - mu[, m]
                y <- sum(x2/sigma[, m] * x2)
                x2 <- x2 %*% beta2
                y <- y - sum((x2 %*% M[[m]]) * x2)
                
                return(y)
            })
            
            loglik[, m] <- (-loglik0[, m] - dt - G * log(2 * PI))/2
        }
        
        loglik <- t(t(loglik) + log(pi))
        sconst <- rowMaxs(loglik)
        loglik2 <- loglik - sconst
        lik <- exp(loglik2)
        
        # compute total loglik
        tot[it] <- sum(log(rowSums(lik)) + sconst)
        
        if (it >= est_z | is.null(z)) {
            z <- lik/rowSums(lik)  # n * M0
        }
        
        nz <- colSums(z)
        
        # get expectation of W
        for (m in 1:M0) {
            Wm[[m]] <- t(M[[m]] %*% t(beta) %*% ((Y - mu[, m])/sigma[, m]))  # n * K0
        }
        
        
        # imputation
        if (it >= impt_it + 1) {
            for (iii in 1:num_mc) {
                if (M0 > 1) {
                  im <- apply(z, 1, function(x) which(rmultinom(1, 1, x) == 1))  # sample membership
                } else {
                  im <- rep(1, n)
                }
                for (i in 1:n) {
                  m <- im[i]
                  vr <- M[[m]]
                  f_i <- rmvnorm(1, Wm[[m]][i, ], vr)
                  
                  ind <- which(Y0[, i] <= cutoff)
                  
                  ms <- (mu[ind, m] + beta[ind, ] %*% t(f_i)) * gene_sd[ind] + gene_mean[ind]
                  sds <- sqrt(sigma[ind, m]) * gene_sd[ind]
                  p <- pg[ind, celltype[i]]  # need celltype
                  
                  prob <- pnorm(cutoff, mean = ms, sd = sds)  # compute x<0 prob
                  prob_drop <- (1 - p)/(prob * p + (1 - p))
                  I_drop <- rbinom(length(ind), 1, prob_drop)
                  
                  
                  # imputation for dropout
                  impt <- rep(0, length(ind))
                  impt[I_drop == 1] <- rnorm(sum(I_drop == 1), ms[I_drop == 1], sds[I_drop == 
                    1])
                  
                  # imputation for non-dropout
                  if (sum(I_drop == 0) > 0) {
                    impt[I_drop == 0] <- rtnorm(sum(I_drop == 0), upper = cutoff, 
                      mean = ms[I_drop == 0], sd = sds[I_drop == 0])
                  }
                  
                  impt[impt > upper] <- upper
                  impt[impt < lower] <- lower
                  
                  Y[ind, i] <- (impt - gene_mean[ind])/gene_sd[ind]
                }
                
                # get expectation of Z
                for (m in 1:M0) {
                  loglik0[, m] <- apply(Y, 2, function(x) {
                    x2 <- x - mu[, m]
                    y <- sum(x2/sigma[, m] * x2)
                    x2 <- x2 %*% beta2
                    y <- y - sum((x2 %*% M[[m]]) * x2)
                    
                    return(y)
                  })
                  
                  loglik[, m] <- (-loglik0[, m] - dt - G * log(2 * PI))/2
                }
                
                loglik <- t(t(loglik) + log(pi))
                sconst <- rowMaxs(loglik)
                loglik2 <- loglik - sconst
                lik <- exp(loglik2)
                
                if (it >= est_z | is.null(z)) {
                  z <- lik/rowSums(lik)  # n * M0
                }
                
                nz <- colSums(z)
                
                # get expectation of W
                for (m in 1:M0) {
                  Wm[[m]] <- t(M[[m]] %*% t(beta) %*% ((Y - mu[, m])/sigma[, m]))  # n * K0
                }
            }
        }
        
        
        Vm <- lapply(1:M0, function(m) M[[m]] * nz[m])
        
        
        # M step
        res <- foreach(g = 1:G, .combine = rbind) %dopar% {
            V <- 0
            for (m in 1:M0) V <- V + Vm[[m]]/sigma[g, m]
            W_temp <- c()
            Y_temp <- c()
            
            for (m in 1:M0) {
                Y_temp <- c(Y_temp, (Y[g, ] - mu[g, m]) * sqrt(z[, m])/sqrt(sigma[g, 
                  m]))
                W_temp <- rbind(W_temp, Wm[[m]] * sqrt(z[, m])/sqrt(sigma[g, m]))
            }
            
            ML <- chol(V)
            
            W_aug <- rbind(W_temp, ML)  # (n+K) * K
            Y_aug <- c(Y_temp, rep(0, K0))  # G*(n+K)
            
            
            penalty <- penl  # sigma^2
            fit1m <- glmnet(W_aug, Y_aug, family = "gaussian", alpha = 1, intercept = F, 
                standardize = F, nlambda = 1, lambda = mean(penalty)/(M0 * n + K0)/var(Y_aug), 
                penalty.factor = rep(1, K0))  # K dimensional, n+K data
            
            nb <- fit1m$beta[, 1]
            
            sg <- sapply(1:M0, function(m) {
                (sum((Y[g, ] - mu[g, m] - nb %*% t(Wm[[m]]))^2 * z[, m]) + sum((nb %*% 
                  Vm[[m]]) * nb))
            })
            c(nb, rep((sum(sg) + 1)/(n + 3), M0))
        }
        
        beta <- matrix(res[, 1:K0], ncol = K0)
        sigma <- matrix(res[, -1:-K0], ncol = M0)
        sigma[sigma > 9] <- 9
        sigma[sigma < 1e-04] <- 1e-04
        
        # max lambda
        if (max_lambda & it >= est_lam) {
            ind_m <- which((nz > K0 + 1) & (nz > 3))  # at least K0 samples to estimate lambda
            mt <- length(ind_m)
            if (mt > 1) {
                Em <- sapply(ind_m, function(m) {
                  # 1:M0
                  colSums(Wm[[m]]^2 * z[, m]) + diag(M[[m]]) * nz[m]
                })
                # print(round(Em,4))
                
                lambda[-ind_m, ] <- 1/M0
                lam_init <- t(Em + 1)/(nz[ind_m] + 1)
                lam_init <- t(t(lam_init)/colSums(lam_init)) * mt/M0
                
                for (k in 1:K0) {
                  if (all(beta[, k] == 0)) {
                    lambda[, k] <- 1/M0
                  } else {
                    lambda[ind_m, k] <- solnp(lam_init[, k], fun = fn, eqfun = equal, 
                      eqB = mt/M0, LB = rep(0, mt), UB = rep(mt/M0, mt), control = list(inner.iter = 500, 
                        trace = 0), nz = nz[ind_m], Em = Em[k, ])$pars
                  }
                }
            } else {
                lambda <- matrix(1/M0, M0, K0)
            }
        }
        
        # maximizme mu
        for (m in 1:M0) {
            sigma_temp <- (nz[m] + sigma[, m]/sigma0)
            beta2 <- beta/sigma_temp
            
            bigM <- diag(1/sigma_temp) - beta2 %*% solve(diag(sigma0/lambda[m, ]) + 
                t(beta2) %*% beta) %*% t(beta2)
            
            mu[, m] <- bigM %*% (Y %*% z[, m])
        }
        
        # max pi
        pi <- (nz + pi_alpha - 1)/(n + pi_alpha * M0 - M0)
        # plot beta every 10 iter == iter
        if (verbose & it%%10 == 0) {
            print(paste(it, tot[it]))
            print("number of cells for each cluster: ")
            print(round(nz))
            print("sum of factor variances for each cluster: ")
            print(round(rowSums(lambda), 2))
            
            if (!is.null(celltype_true)) {
                print("clustering randInd: ")
                print(mclust::adjustedRandIndex(im, celltype_true))
            }
            
        }
    }
    
    Y <- Y * gene_sd + gene_mean
    
    return(list(loglik = tot, pi = pi, mu = mu, sigma = sigma, beta = beta, lambda = lambda, 
        z = z, Ef = Wm, Varf = M, Y = Y, geneM = gene_mean, geneSd = gene_sd))
}
