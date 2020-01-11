#' Sampling multiple copies of imputed values and cell factors at given parameters.
#'
#' \code{do_impute} get the posterior mean and variance for the imputed values and factors given the parameters, which can be returned from \emph{SIMPLE} or \emph{SIMPLE_B}.
#'
#' @param dat scRNASeq data matrix. Each row is a gene, each column is a cell.
#' @param Y  Initial imputed data, output from \emph{SIMPLE_B} or \emph{SIMPLE}.
#' @param beta  Factor loadings, output from \emph{SIMPLE_B} or \emph{SIMPLE}.
#' @param lambda  Variances of factors for each cluster, output from \emph{SIMPLE_B} or \emph{SIMPLE}.
#' @param sigma  Variances of idiosyncratic noises for each cluster, output from \emph{SIMPLE_B} or \emph{SIMPLE}.
#' @param mu  Mean expression for each cluster, output from \emph{SIMPLE_B} or \emph{SIMPLE}.
#' @param pi  Probabilites of cells belong to each cluster, output from \emph{SIMPLE_B} or \emph{SIMPLE}.
#' @param pos_mean  Gene mean. If centerized each gene before estimating the parameters, provide the overall mean of gene expression removed from the data matrix, output from \emph{SIMPLE_B} or \emph{SIMPLE}.
#' Default is NULL.
#' @param pos_sd  Gene standard deviation. If scaled each gene before estimating the parameters, provide the overall standard deviation of gene expression removed from the data matrix, output from \emph{SIMPLE_B} or \emph{SIMPLE}. Default is NULL.
#' Default is NULL.
#'
#' @return \code{do_impute} returns a list of imputation results in the following order.
#' \enumerate{
#'   \item{loglik}{The log-likelihood of the imputed gene expression at each iteration.}
#'   \item{impt}{A matrix contains the expectation of imputed expression.}
#'  \item{impt_var}{A matrix contains the variance of each imputed entry.}
#'   \item{EF}{Posterior means of factors}
#'   \item{varF}{Posterior covariance matrix of factors}
#'  \item{consensus_cluster}{Score for the clustering stability of each cell by multiple imputations. }
#' }
#' @author Zhirui Hu, \email{zhiruihu@g.harvard.edu}
#' @author Songpeng Zu, \email{songpengzu@g.harvard.edu}
#'
do_impute <- function(dat, Y, beta, lambda, sigma, mu, pi, pos_mean = NULL, pos_sd = NULL, 
    celltype = NULL, mcmc = 10, burnin = 2, verbose = F, pg = 0.5, cutoff = 0.1) {
    PI = 3.14159265359
    # initiation
    G <- nrow(Y)
    n <- ncol(Y)
    M0 <- ncol(mu)
    K0 <- ncol(beta)
    
    if (is.null(pos_mean)) 
        pos_mean <- rep(0, G)
    if (is.null(pos_sd)) 
        pos_sd <- rep(1, G)
    
    if (is.null(celltype)) 
        celltype <- rep(1, n)
    MB <- max(celltype)
    
    if (length(pg) == 1) 
        pg <- matrix(pg, G, MB)
    
    Y <- (Y - pos_mean)/pos_sd
    
    # A = diag(1,K0,K0) store result
    loglik <- matrix(0, n, M0)
    tot <- rep(0, mcmc + burnin)
    record_impt <- matrix(0, G, n)
    record_impt2 <- matrix(0, G, n)
    record_EF <- matrix(0, n, K0)
    record_F2 <- matrix(0, n, K0)
    record_varF <- matrix(0, n, K0^2)
    
    consensus_cluster <- matrix(0, n, n)
    
    M <- list()  # (BSigma^-1B + Lambda^-1)^-1
    Wm <- list()  # mean of each factor, K0*n
    
    
    for (it in 1:(mcmc + burnin)) {
        # sample Z
        for (m in 1:M0) {
            beta2 <- beta/sigma[, m]
            M[[m]] <- solve(t(beta) %*% beta2 + diag(1/lambda[m, ]))
            dt <- sum(log(sigma[, m])) + sum(log(lambda[m, ])) - log(det(M[[m]]))  # !!! logdet !=det(log=T)
            loglik[, m] <- apply(Y, 2, function(x) {
                x2 <- x - mu[, m]
                y <- sum(x2/sigma[, m] * x2)
                x2 <- x2 %*% beta2
                y <- y - sum((x2 %*% M[[m]]) * x2)
                
                return((-y - dt - G * log(2 * PI))/2)
            })
        }
        
        loglik <- t(t(loglik) + log(pi))
        sconst <- rowMaxs(loglik)
        loglik2 <- loglik - sconst
        lik <- exp(loglik2)
        
        # compute total loglik
        tot[it] <- sum(log(rowSums(lik)) + sconst)
        
        z <- lik/rowSums(lik)  # n * M0
        nz <- colSums(z)
        
        ## overall consensus clustering
        if (it > burnin) 
            consensus_cluster <- consensus_cluster + z %*% t(z)
        
        # print(nz) print(tot[it])
        
        
        # get expectation of W
        for (m in 1:M0) {
            Wm[[m]] <- t(M[[m]] %*% t(beta) %*% ((Y - mu[, m])/sigma[, m]))  # n * K0
        }
        
        # imputation
        im <- apply(z, 1, function(x) which(rmultinom(1, 1, x) == 1))  # sample membership
        # Y = foreach(i=1:n, .combine = cbind) %dopar% {
        for (i in 1:n) {
            m <- im[i]
            vr <- M[[m]]
            f_i <- rmvnorm(1, Wm[[m]][i, ], vr)
            
            ind <- which(dat[, i] <= cutoff)
            
            ms <- (mu[ind, m] + beta[ind, ] %*% t(f_i)) * pos_sd[ind] + pos_mean[ind]
            sds <- sqrt(sigma[ind, m]) * pos_sd[ind]
            p <- pg[ind, celltype[i]]
            prob <- pnorm(cutoff, mean = ms, sd = sds)  # compute x<0 prob
            prob_drop <- (1 - p)/(prob * p + (1 - p))
            I_drop <- rbinom(length(ind), 1, prob_drop)
            
            # imputation for dropout
            impt <- rep(0, length(ind))
            impt[I_drop == 1] <- rnorm(sum(I_drop == 1), ms[I_drop == 1], sds[I_drop == 
                1])
            
            # imputation for non-dropout
            if (sum(I_drop == 0) > 0) {
                impt[I_drop == 0] <- rtnorm(sum(I_drop == 0), upper = cutoff, mean = ms[I_drop == 
                  0], sd = sds[I_drop == 0])
            }
            
            Y[ind, i] <- (impt - pos_mean[ind])/pos_sd[ind]
            # return(res)
        }
        
        # record
        if (it > burnin) {
            record_impt <- record_impt + Y
            record_impt2 <- record_impt2 + Y^2
            
            for (m in 1:M0) {
                record_EF <- record_EF + Wm[[m]] * z[, m]
                record_F2 <- record_F2 + Wm[[m]]^2 * z[, m]
                record_varF <- record_varF + z[, m] %*% t(c(M[[m]]))
            }
        }
        # print loglik
        if (verbose & it%%20 == 0) {
            print(paste(it, tot[it]))
            print(nz)
        }
    }
    
    varF <- record_F2/mcmc - (record_EF/mcmc)^2 + record_varF[, seq(1, K0^2, by = (K0 + 
        1))]/mcmc
    return(list(loglik = tot, impt = (record_impt/mcmc) * pos_sd + pos_mean, impt_var = record_impt2/mcmc - 
        (record_impt/mcmc)^2, EF = record_EF/mcmc, varF = varF, consensus_cluster = consensus_cluster/mcmc))
}
