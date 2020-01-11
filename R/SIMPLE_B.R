#' Initialize imputation for each individual gene with known dropout rate for each cell type
#'
#' @details
#' Fit each gene expression by a zero-inflated censored Gaussian distribution and return a random sample of imputed values for initialization
#' @param Y2  scRNASeq data matrix. Each row is a gene, each column is a cell.
#' @param clus A numeric vector for the cell type labels of cells in the scRNASeq. The labels must start from 1 to the number of types.
#' @param bulk Gene expression matrix of bulk RNASeq, only used for plot if verbose = T.
#' @param pg1 lower bound of Dropout rate (at least 1- pg1) for each cell type.
#' @param cutoff The value below cutoff is treated as no expression. Default = 0.1.
#' @param verbose Whether to plot some intermediate result. Default = False.
#'
#' @return Imputed gene expression matrix, treat each gene independently
#' @import MASS
#' @import glmnet
#' @import matrixStats
#' @import Rsolnp
#' @importFrom mixtools rmvnorm
#' @importFrom msm rtnorm
#'
#' @author Zhirui Hu, \email{zhiruihu@g.harvard.edu}
#' @author Songpeng Zu, \email{songpengzu@g.harvard.edu}

init_impute_bulk <- function(Y2, clus, bulk, pg1, cutoff = 0.1, verbose = F) {
    impute0 <- Y2
    G <- nrow(Y2)
    pg <- pg1  # matrix(0,G, M0)
    # ndropout = pg1
    M0 <- max(clus)
    for (i in 1:M0) {
        temp_dat <- as.matrix(Y2[, clus == i])
        # p0 = rowMeans(temp_dat <= cutoff)
        impute <- temp_dat
        
        result <- ztruncnorm(temp_dat, cutoff = cutoff, p_min = 0.01, p_max = pg1[, 
            i])  # pg1 >= 0.1
        mu <- result[[1]][, 1]
        sd <- result[[1]][, 2]
        pg[, i] <- result[[2]]
        
        
        ind <- which(is.na(mu))
        if (verbose) 
            print(paste(length(ind), "not fit TN"))
        mu[ind] <- rowMeans(Y2[ind, clus == i, drop = F])
        sd[ind] <- rowSds(Y2[ind, clus == i, drop = F])
        pg[ind, i] <- pg1[ind, i]  # 1
        
        
        if (verbose) {
            par(mfrow = c(2, 2), pch = 16)
            plot(bulk[, i], mu, col = rgb(1, 0, 0, 0.5), xlab = "bulk", ylab = "est_mu", 
                main = paste("cluster", i))
            a <- coef(lm(mu ~ bulk[, i], weights = bulk[, i]^2))
            abline(a, col = 2)
            
            plot(bulk[, i], rowMeans(temp_dat), col = rgb(0, 1, 0, 0.5), xlab = "bulk", 
                ylab = "raw_mean")
            b <- coef(lm(rowMeans(temp_dat) ~ bulk[, i], weights = bulk[, i]^2))
            abline(b, col = 2)
            
            # dropout rate
            plot(pg1[, i], pg[, i], col = rgb(0, 1, 1, 0.5), xlab = "amplified rate", 
                ylab = "prob normal component")
        }
        
        
        # imputation for each gene
        for (j in 1:G) {
            imp <- which(temp_dat[j, ] <= cutoff)
            if (length(imp) == 0) 
                next
            prob <- pnorm(cutoff, mean = mu[j], sd = sd[j])  # compute x<0 prob
            
            # imputation for dropout from norm, zero comp, censored part q * d, pg1 = 1-d, p
            # = (1-d)*q, q = pg/pg1
            prob_drop <- (pg[j, i]/pg1[j, i] - pg[j, i])/(prob * pg[j, i] + (1 - 
                pg[j, i]))
            # 1-q
            prob_zero <- (1 - pg[j, i]/pg1[j, i])/(prob * pg[j, i] + (1 - pg[j, i]))
            # q(1-d)*prob
            prob_trunc <- prob * pg[j, i]/(prob * pg[j, i] + (1 - pg[j, i]))
            
            I_drop <- rmultinom(length(imp), 1, c(prob_drop, prob_zero, prob_trunc))
            
            # imputation for dropout
            impute[j, imp[I_drop[1, ] == 1]] <- rnorm(sum(I_drop[1, ] == 1), mu[j], 
                sd[j])
            
            # imputation for others
            impute[j, imp[I_drop[2, ] == 1]] <- rnorm(sum(I_drop[2, ] == 1), -1, 
                0.5)
            
            impute[j, imp[I_drop[3, ] == 1]] <- rtnorm(sum(I_drop[3, ] == 1), upper = cutoff, 
                mean = mu[j], sd = sd[j])
        }
        
        # print(impute[j,imp])
        
        
        impute0[, clus == i] <- impute
        
        if (verbose) {
            plot(bulk[, i], rowMeans(impute), col = rgb(0, 0, 1, 0.5), xlab = "bulk", 
                ylab = "imputed mean")
            abline(b, col = 2)
        }
    }
    
    impute0[is.na(impute0)] <- 0
    
    
    return(impute0)
    
}

#' Impute zero entries and clustering for scRNASeq data integrating bulk RNASeq
#'
#' Impute zero entries in the gene expression matrix based on the expression level in the similar cells and gene-gene correlation.
#' Zero entries in the observed expression matrix come from either molecule loss during the experiment ('dropout') or too low expression to be measured.
#' We use Monte Carlo EM algorithm to sample the imputed values and learn the parameters.
#' We use this function to further integrate the bulk RNASeq data for the similar cell types as the scRNASeq to infer the dropout rate.
#'
#' @details
#' We assume that the cells come from M0 clusters. Within each cell cluster, the
#'   'true' gene expression is modeled by a multivariate Gaussian distribution
#'   whose covariance matrix can be composed into a low rank matrix (a couple of
#'   latent gene modules) and idiosyncratic noises. Gene modules are shared
#'   among cell clusters though the coexpression level of each gene module can
#'   be different. Suppose there are G genes and n cells. For each cell cluster,
#'   the gene expression follows \eqn{Y|Z=m~MVN(\mu_m, B\Lambda_m B^T +
#'   \Sigma_m)} where B is a G by K0 matrix, \eqn{\Sigma_m} is a G by G diagonal
#'   matrix whose diagonal entries are specified by \emph{sigma}, and
#'   \eqn{\Lambda_m} is a K0 by K0 diagonal matrix whose diagonal entries are
#'   specified by \emph{lambda}. \eqn{P(Z_m) = \pi_m} where
#'   \eqn{\pi~Dir(\alpha)}. 
#'   
#'   The algorithm first runs Monte Carlo EM using only
#'   the genes with low dropout rate (initial phase) and initializes factor
#'   loadings and clustering membership. Then it runs another rounds of Monte
#'   Carlo EM using all the genes. In the initial phase, we use the genes with dropout rate less than
#' \emph{1 - p_min}; if the number of genes is less than \emph{min_gene}, we
#' rank the genes by the number cells with nonzero expression and keep the top
#' \emph{min_gene} genes. If \emph{fix_num} is true, then we always keep the top 
#' \emph{min_gene} genes in the initial phase.

#' @param dat scRNASeq data matrix. Each row is a gene, each column is a cell.
#' @param K0 Number of latent gene modules. See details.
#' @param bulk Bulk RNASeq data matrix. Each row is a gene which must be ordered the same as the scRNASeq data. Each column is a cell type which must be ordered as the cell type label in \emph{celltype}.
#' @param celltype A numeric vector for labels of cells in the scRNASeq. Each cell type has corresponding mean expression in the bulk RNASeq data. The labels must start from 1 to the number of types. If NULL, all cells are treated as a single cell type and the input bulk RNASeq should also have one column that is the mean expression over all the cell types.
#' @param M0 Number of clusters. See details.
#' @param clus Initial clustering of scRNASeq data. If NULL, the function will use PCA and Kmeans to do clustering initially.
#' @param K The number of PCs used in the initial clustering. Default = 20.
#' @param iter The number of EM iterations using full data set. See details.
#' @param est_z The iteration starts to update Z.
#' @param impt_it The iteration starts to sample new imputed values in initial phase. See details.
#' @param max_lambda Whether to maximize over lambda.
#' @param est_lam The iteration starts to estimate lambda.
#' @param penl L1 penalty for the factor loadings.
#' @param sigma0 The variance of the prior distribution of \eqn{\mu}.
#' @param pi_alpha The hyperparameter of the prior distribution of \eqn{\pi}. See details.
#' @param beta A G by K0 matrix. Initial values for factor loadings (B). If null, beta will initialze from normal distribution with mean zero and variance M0/K0. See details.
#' @param lambda A M0 by K0 matrix. Initial values for the variances of factors. Each column is for a cell cluster. If null, lambda will initialize to be 1/M0. See details.
#' @param sigma A G by M0 matrix. Initial values for the variance of idiosyncratic noises. Each column is for a cell cluster. If null, sigma will initialize to be 1. See details.
#' @param mu A G by M0 matrix. Initial values for the gene expression mean of each cluster. Each column is for a cell cluster. If NULL, it will take the sample mean of cells weighted by the probability in each cluster. See details.
#' @param p_min Initialize parameters using genes expressed in at least \emph{p_min} proportion of cells. If the number genes selected is less than \emph{min_gene}, select \emph{min_gene} genes with higest proportion of non zeros. Default = 0.8.
#' @param min_gene Initialize parameters using genes expressed in at least \emph{p_min} proportion of cells. If the number genes selected is less than \emph{min_gene}, select \emph{min_gene} genes with higest proportion of non zeros. Default = 1000.
#' @param fix_num If true, always use top \emph{min_gene} genes with higest proportion of non zeros in the initial phase. Default  = F. See details.
#' @param cutoff The value below cutoff is treated as no expression. Default = 0.1.
#' @param verbose Whether to show some intermediate results. Default = False.
#' @param num_mc The number of Gibbs steps for generating imputed data when the parameters are updated during Monte Carlo EM. Default = 3.
#' @param mcmc The number of Gibbs steps to sample imputed data after EM. Default = 50. 
#' @param burnin The number of burnin steps before sample imputed data after EM. Default = 2.
#' @return \code{SIMPLE_B} returns a list of results in the following order.
#' \enumerate{
#' \item{loglik}{The log-likelihood of the imputed gene expression at each iteration.}
#' \item{pi}{Probabilites of cells belong to each cluster.}
#' \item{mu}{Mean expression for each cluster}
#' \item{sigma}{Variances of idiosyncratic noises for each cluster.}
#' \item{beta}{Factor loadings.}
#' \item{lambda}{Variances of factors for each cluster.}
#' \item{z}{The probability of each cell belonging to each cluster.}
#' \item{Yimp0}{A matrix contains the expectation of imputed expression.}
#' \item{pg}{A G by M0 matrix, dropout rate for each gene in each cluster.}
#'  \item{impt}{A matrix contains the mean of each imputed entry by sampling multiple imputed values at MLE. If mcmc <= 0, output imputed expressoin matrix at last step of EM}
#'  \item{impt_var}{A matrix contains the variance of each imputed entry by sampling multiple imputed values at MLE. NULL if mcmc <= 0.}
#'     \item{Ef} {If mcmc >0, output posterior means of factors
#'     given observed data (a n by K0 matrix). If mcmc <= 0, output conditional expectation of the factors for each cluster \eqn{E(f_i|z_i= m)} 
#'    at the last step of EM. A list with length M0, 
#'    each element in the list is a n by K0 matrix.}
#'     \item{Varf} {If mcmc >0, output posterior variances of
#'     factors given observed data (a n by K0 matrix). If mcmc <= 0, output conditional covariance matrix of factors for each cluster \eqn{Var(f_i|z_i = m)} at the last step of EM. 
#'      A list with length M0, each element in the list is a K0 by K0 matrix.}
#'  \item{consensus_cluster}{Score for the clustering stability of each cell by multiple imputations. NULL if mcmc <=0 }
#' }
#' @import doParallel
#' @importFrom foreach foreach
#' @seealso \code{\link{SIMPLE}}
#' @examples
#' library(foreach) 
#' library(doParallel) 
#' library(SIMPLE) 
#' 
#' M0 = 3  # simulate number of clusters 
#' n = 300 # number of cells 
#' 
#' simu_data = simulation_bulk(n=300, S0 = 20, K = 6, MC=M0, block_size = 32, indepG = 1000 - 32*6, verbose=F, overlap=0) 
#' 
#' Y2 = simu_data$Y2 
#' K0 = 6 # number of factors 
#' 
#' registerDoParallel(cores = 6)  # parallel 
#' 
#' # estimate the parameters, only input the overall mean of all cell types from bulk 
#' result <- SIMPLE_B(Y2, K0, data.frame(simu_data$bulk), M0, celltype=rep(1, n), clus = NULL, K = 20, p_min = 0.5, min_gene = 200,cutoff=0.01, max_lambda=T) 
#' 
#' # evaluate cluster performance 
#' celltype_true = simu_data$Z 
#' mclust::adjustedRandIndex(apply(result$z,1, which.max), celltype_true) 
#' # or redo clustering based on imputed values (sometimes work better for real data) 
#' getCluster(result$impt, celltype_true, Ks = 20, M0 = M0)[[1]] 
#'
#' @author Zhirui Hu, \email{zhiruihu@g.harvard.edu}
#' @author Songpeng Zu, \email{songpengzu@g.harvard.edu}
#' @export

SIMPLE_B <- function(dat, K0, bulk, celltype, M0 = 1, clus = NULL, K = 20, 
    iter = 10, est_z = 1, impt_it = 3, max_lambda = F, est_lam = 1, penl = 1, sigma0 = 100, 
    pi_alpha = 1, beta = NULL, lambda = NULL, sigma = NULL, mu = NULL, p_min = 0.8, 
    min_gene = 300, cutoff = 0.1, verbose = F, num_mc = 3, fix_num = F, mcmc = 50, 
    burnin = 2) {
    # EM algorithm initiation
    G <- nrow(dat)
    n <- ncol(dat)
    Y <- dat  # imputed matrix
    # remove gene mean for factor, otherwise may be confounded with B, i.e. u = B1,
    # f1=rep(1,n)
    gene_mean <- rep(0, G)
    
    # prob of z
    pi <- rep(1/M0, M0)
    
    # random start sum M0
    if (is.null(lambda)) 
        lambda <- matrix(1/M0, M0, K0)
    
    if (is.null(mu)) 
        mu <- matrix(0, G, M0)
    if (is.null(sigma)) 
        sigma <- matrix(1, G, M0)
    
    if (is.null(beta)) 
        beta <- matrix(rnorm(G * K0), G, K0)/sqrt(K0) * sqrt(M0)
    
    
    # inital impution only for low dropout genes
    n1 <- rowMeans(dat > cutoff)
    if (fix_num) {
        hq_ind <- order(n1, decreasing = T)[1:min_gene]
    } else {
        hq_ind <- which(n1 >= p_min)
        if (length(hq_ind) < min_gene) 
            hq_ind <- order(n1, decreasing = T)[1:min_gene]  # need to change back
    }
    # regression for hq genes and estimate dropout rate for all genes
    MB <- max(celltype)
    pg <- matrix(0, G, MB)
    for (i in 1:MB) {
        # b = coef(lm(rowMeans(dat[hq_ind, celltype==i])~bulk[hq_ind,i], weights =
        # bulk[hq_ind,i]^2))
        b <- c(0, 1)
        pg[, i] <- (rowMeans(dat[, celltype == i]) + p_min)/(1 + b[1] + b[2] * bulk[, 
            i])
    }
    
    pg[pg > 1] <- 1
    pg[pg < 0.1] <- 0.1
    
    message(paste("initial impution for hq genes: ", length(hq_ind)))  # low dropout
    
    res <- init_impute_bulk(dat[hq_ind, ], celltype, bulk[hq_ind, , drop = F], pg[hq_ind, 
        , drop = F], cutoff = cutoff, verbose = F)  # verbose
    
    # init clustering
    if (is.null(clus)) {
        
        Y2_scale <- t(scale(t(res)))
        s <- svd(Y2_scale)
        # if (clus_opt == 1) {
        km0 <- kmeans(t(Y2_scale) %*% s$u[, 1:K], M0, iter.max = 80, nstart = 300)
        # } else { km0 <- kmeans(s$v[, 1:K], M0, iter.max = 80, nstart = 300) }
        clus <- km0$cluster
        
        if (verbose) {
            print(xtabs(~clus))
            print("initial clustering randInd: ")
            
            if (!is.null(celltype)) {
                print(mclust::adjustedRandIndex(clus, celltype))  #_true
            }
        }
    }
    
    z <- matrix(0, n, M0)
    for (m in 1:M0) z[clus == m, m] <- 1
    
    message("initial estimate factors: ")
    impute_hq <- EM_impute(res, dat[hq_ind, ], pg[hq_ind, , drop = F], M0, K0, cutoff, 
        30, beta[hq_ind, ], sigma[hq_ind, , drop = F], lambda, pi, z, mu = NULL, 
        celltype = celltype, penl, est_z, max_lambda, est_lam, impt_it, sigma0, pi_alpha, 
        verbose = verbose, num_mc = num_mc)  # iter = 30, mu = NULL
    
    
    beta[hq_ind, ] <- impute_hq$beta
    sigma[hq_ind, ] <- impute_hq$sigma
    mu[hq_ind, ] <- impute_hq$mu  # check this, as imputed mu may not be zero
    Y[hq_ind, ] <- impute_hq$Y
    gene_mean[hq_ind] <- impute_hq$geneM
    z <- impute_hq$z
    
    
    nz <- colSums(impute_hq$z)
    Vm <- lapply(1:M0, function(m) impute_hq$Varf[[m]] * nz[m])
    
    # inital beta for other genes
    message("initial estimate beta for lq genes:")
    
    lq_ind <- setdiff(1:G, hq_ind)
    # estimate beta and impute: only for positive part? (only impute for genes with
    # more than 10% nonzero)
    
    # M step also estimate mu sapply(lq_ind, function(g){#
    res <- foreach(g = lq_ind, .combine = rbind) %dopar% {
        V <- 0
        for (m in 1:M0) V <- V + Vm[[m]]/sigma[g, m]
        W_temp <- c()
        Y_temp <- c()
        
        for (m in 1:M0) {
            Y_temp <- c(Y_temp, dat[g, ] * sqrt(z[, m])/sqrt(sigma[g, m]))
            Wb <- impute_hq$Ef[[m]] * sqrt(z[, m])/sqrt(sigma[g, m])
            Wmu <- matrix(0, n, M0)
            Wmu[, m] <- sqrt(z[, m])/sqrt(sigma[g, m])  # n * M
            W_temp <- rbind(W_temp, cbind(Wmu, Wb))
        }
        
        ML <- cbind(matrix(0, K0, M0), chol(V))
        
        W_aug <- rbind(W_temp, ML)  #(n+K) * (M + K)
        Y_aug <- c(Y_temp, rep(0, K0))  #G*(n+K)
        
        penalty <- penl/(M0 * n + K0)/var(Y_aug)  # sigma^2
        fit1m = glmnet(W_aug, Y_aug, family = "gaussian", alpha = 1, intercept = F, 
            standardize = F, nlambda = 1, lambda = penalty * K0/(M0 + K0), penalty.factor = c(rep(0, 
                M0), rep(1, K0)))  #K dimensional, n+K data
        
        coeff = fit1m$beta[-1:-M0, 1]
        tempmu = fit1m$beta[1:M0, 1]
        
        sg = sapply(1:M0, function(m) {
            (sum((dat[g, ] - tempmu[m] - coeff %*% t(impute_hq$Ef[[m]]))^2 * z[, 
                m]) + sum((coeff %*% Vm[[m]]) * coeff))
        })
        c(fit1m$beta[, 1], rep((sum(sg) + 1)/(n + 3), M0))
    }
    
    mu[lq_ind, ] = matrix(res[, 1:M0], ncol = M0)
    beta[lq_ind, ] = matrix(res[, (M0 + 1):(K0 + M0)], ncol = K0)
    
    sigma[lq_ind, ] = matrix(res[, -1:-(K0 + M0)], ncol = M0)
    sigma[sigma > 9] = 9
    sigma[sigma < 1e-04] = 1e-04
    
    
    
    # imputation set dropout rate as pg
    if (M0 > 1) {
        im = apply(z, 1, function(x) which(rmultinom(1, 1, x) == 1))  # sample membership
    } else {
        im = rep(1, n)
    }
    for (i in 1:n) {
        m = im[i]
        
        ind = which(dat[lq_ind, i] <= cutoff)
        ind = lq_ind[ind]
        
        ms = mu[ind, m] + beta[ind, , drop = F] %*% impute_hq$Ef[[m]][i, ]
        sds = sqrt(sigma[ind, m])
        p = pg[ind, celltype[i]]  # need celltype
        
        prob = pnorm(cutoff, mean = ms, sd = sds)  # compute x<0 prob
        prob_drop = (1 - p)/(prob * p + (1 - p))
        I_drop = rbinom(length(ind), 1, prob_drop)
        
        
        # imputation for dropout
        impt = rep(0, length(ind))
        impt[I_drop == 1] = rnorm(sum(I_drop == 1), ms[I_drop == 1], sds[I_drop == 
            1])
        
        # imputation for non-dropout
        if (sum(I_drop == 0) > 0) {
            impt[I_drop == 0] = rtnorm(sum(I_drop == 0), upper = cutoff, mean = ms[I_drop == 
                0], sd = sds[I_drop == 0])
        }
        
        Y[ind, i] = impt
    }
    
    impute <- matrix(0, n, G)
    
    # EM for all genes
    message("impute for all genes")
    impute_result <- EM_impute(Y, dat, pg, M0, K0, cutoff, iter, beta, sigma, impute_hq$lambda, 
        impute_hq$pi, impute_hq$z, mu = NULL, celltype = celltype, penl, est_z, max_lambda, 
        est_lam, impt_it = 1, sigma0, pi_alpha, verbose = verbose, num_mc = num_mc)
    
    for (m in 1:M0) {
        impute <- impute + t(impute_result$mu[, m] + impute_result$beta %*% t(impute_result$Ef[[m]])) * 
            impute_result$z[, m]
    }
    
    impute <- t(impute) * impute_result$geneSd + impute_result$geneM
    
    if (mcmc > 0) {
        message("multiple impution sampling")
        result2 = do_impute(dat, impute_result$Y, impute_result$beta, impute_result$lambda, 
            impute_result$sigma, impute_result$mu, impute_result$pi, impute_result$geneM, 
            impute_result$geneSd,celltype, mcmc = mcmc, burnin = burnin, pg = pg, cutoff = cutoff)

        return(list(loglik = impute_result$loglik, pi = impute_result$pi, mu = impute_result$mu, 
            sigma = impute_result$sigma, beta = impute_result$beta, lambda = impute_result$lambda, 
            z = impute_result$z, Yimp0 = impute, pg = pg, impt = result2$impt, impt_var = result2$impt_var, 
            Ef = result2$EF, Varf = result2$varF, consensus_cluster = result2$consensus_cluster))
    } else {
        return(list(loglik = impute_result$loglik, pi = impute_result$pi, mu = impute_result$mu, 
            sigma = impute_result$sigma, beta = impute_result$beta, lambda = impute_result$lambda, 
            z = impute_result$z, Yimp0 = impute, pg = pg, impt = impute_result$Y, 
            impt_var = NULL, Ef = impute_result$Ef, 
            Varf = impute_result$Varf, consensus_cluster = NULL))
    }
}
