#' Initialize imputation for each individual gene
#'
#' @details
#' Fit each gene expression by a zero-inflated censored Gaussian distribution and return a random sample of imputed values for initialization
#' @param Y2  scRNASeq data matrix. Each row is a gene, each column is a cell.
#' @param M0 number of cell types
#' @param clus A numeric vector for the cell type labels of cells in the scRNASeq. The labels must start from 1 to the number of types (M0).
#' @param p_min Restrict the max dropout rate to be 1-p_min. Default = 0.6.
#' @param cutoff The value below cutoff is treated as no expression. Default = 0.1.
#' @param verbose Whether to plot some intermediate result. Default = False.
#'
#' @return Imputed gene expression matrix, treat each gene independently
#' @author Zhirui Hu, \email{zhiruihu@g.harvard.edu}
#' @author Songpeng Zu, \email{songpengzu@g.harvard.edu}
init_impute <- function(Y2, M0, clus, p_min = 0.6, cutoff = 0.1, verbose = F) {
    # fit truncnorm for each cluster
    impute <- Y2
    G <- nrow(Y2)
    pg <- matrix(0, G, M0)
    for (i in 1:M0) {
        temp_dat <- as.matrix(Y2[, clus == i])
        result <- ztruncnorm(temp_dat, cutoff = cutoff, p_min = p_min)
        mu <- result[[1]][, 1]
        sd <- result[[1]][, 2]
        pg[, i] <- result[[2]]
        
        pg[pg > 0.99] <- 0.99
        
        ind <- which(is.na(mu))
        if (verbose) 
            print(paste(length(ind), "not fit TN"))
        mu[ind] <- rowMeans(Y2[ind, clus == i, drop = F])
        sd[ind] <- rowSds(Y2[ind, clus == i, drop = F])
        pg[ind, i] <- 1
        
        if (verbose) 
            print(summary(sd))
        
        for (j in which(clus == i)) {
            ind <- which(Y2[, j] <= cutoff)
            # impute Y2
            ms <- mu[ind]
            sds <- sd[ind]
            p <- pg[ind, i]
            # compute x < - prob
            prob <- pnorm(cutoff, mean = ms, sd = sds)
            prob_drop <- (1 - p)/(prob * p + (1 - p))
            I_drop <- rbinom(length(prob_drop), 1, prob_drop)
            
            # imputation for dropout
            impute[ind[I_drop == 1], j] <- rnorm(sum(I_drop == 1), ms[I_drop == 1], 
                sds[I_drop == 1])
            # imputation for non-dropout
            if (sum(I_drop == 0) > 0) {
                # ms[I_drop==0] - sds[I_drop==0] * dnorm(r)/pnorm(r)
                impute[ind[I_drop == 0], j] <- rtnorm(sum(I_drop == 0), upper = cutoff, 
                  mean = ms[I_drop == 0], sd = sds[I_drop == 0])
            }
        }
    }
    
    impute[is.na(impute)] <- 0
    if (verbose) 
        print(summary(c(impute[Y2 == 0])))
    return(list(impute, pg))
}

#' SIMPLE: Imputing zero entries and clustering for scRNASeq data.
#'
#' \code{SIMPLE} imputes zeros in the gene expression data using the expression level in
#' similar cells and gene-gene correlation. Zero entries in the observed expression matrix
#' come from molecule loss during the experiment ('dropout') or too low expression to
#' be measured. We used Monte Carlo EM algorithm to sample the imputed values and
#' obtain MLEs of all the parameters.
#'
#' @details
#' We assume that the cells come from M0 clusters. Within each cell
#' cluster, the 'true' gene expression is modeled by a multivariate Gaussian
#' distribution whose covariance matrix can be composed into a low rank matrix
#' (a couple of latent gene modules) and idiosyncratic noises. Gene modules are
#' shared among cell clusters, although the coexpression level of each gene module
#' in different cell cluster can be different. \cr
#' Suppose there are G genes and n cells. For each cell
#' cluster, the gene expression follows \eqn{Y|Z=m~MVN(\mu_m, B\Lambda_m B^T +
#' \Sigma_m)} where B is a G by K0 matrix, \eqn{\Sigma_m} is a G by G diagonal
#' matrix whose diagonal entries are specified by \emph{sigma}, and
#' \eqn{\Lambda_m} is a K0 by K0 diagonal matrix whose diagonal entries are
#' specified by \emph{lambda}. \eqn{P(Z_m) = \pi_m} where \eqn{\pi~Dir(\alpha)}.
#'
#' The algorithm first runs Monte Carlo EM using only the genes with low dropout
#' rate (initial phase) and initializes factor loadings and clustering
#' membership. Then it runs more rounds of Monte Carlo EM using all the
#' genes. In the initial phase, we use the genes with dropout rate less than
#' \emph{1 - p_min}; if the number of such genes is less than \emph{min_gene}, we
#' rank the genes by the number cells with nonzero expression and keep the top
#' \emph{min_gene} genes. If \emph{fix_num} is true, then we always keep the top
#' \emph{min_gene} genes in the initial phase.
#'
#' After Monte Carlo EM, we obtain MLE of B, \eqn{\Lambda_m}, etc. If \emph{num_mc} > 0, this function will sample multiple imputed values and cell factors at the MLEs of all the parameters by Gibbs sampling.
#' Based on multiple imputed values, it will evaluate cluster stability for each cell (\emph{consensus_cluster}).
#' It will also output the posterior mean and variance for the imputed values and cell factors.
#'
#' @param dat scRNASeq data matrix. Each row is a gene, each column is a cell.
#' @param K0 Number of latent gene modules. See details.
#' @param M0 Number of clusters. See details.
#' @param clus Initial clustering of scRNASeq data. If NULL, the function will
#'   use PCA and Kmeans to do clustering initially.
#' @param K The number of PCs used in the initial clustering. Default = 20.
#' @param iter The number of EM iterations using full data set. See details.
#' @param est_z The iteration starts to update Z.
#' @param impt_it The iteration starts to sample new imputed values in initial phase. See details.
#' @param max_lambda Whether to maximize over lambda.
#' @param est_lam The iteration starts to estimate lambda.
#' @param penl L1 penalty for the factor loadings.
#' @param sigma0 The variance of the prior distribution of \eqn{\mu}.
#' @param pi_alpha The hyperparameter of the prior distribution of \eqn{\pi}.
#'   See details.
#' @param beta A G by K0 matrix. Initial values for factor loadings (B). If
#'   null, beta will be initialized from normal distribution with mean zero and
#'   variance M0/K0. See details.
#' @param lambda A M0 by K0 matrix. Initial values for the variances of factors.
#'   Each column is for a cell cluster. If null, lambda will initialize to be
#'   1/M0. See details.
#' @param sigma A G by M0 matrix. Initial values for the variance of
#'   idiosyncratic noises. Each column is for a cell cluster. If null, sigma
#'   will initialize to be 1. See details.
#' @param mu A G by M0 matrix. Initial values for the gene expression mean of
#'   each cluster. Each column is for a cell cluster. If NULL, it will take the
#'   sample mean of cells weighted by the probability in each cluster. See
#'   details.
#' @param p_min Initialize parameters using genes expressed in at least
#'   \emph{p_min} proportion of cells. If the number genes selected is less than
#'   \emph{min_gene}, select \emph{min_gene} genes with higest proportion of non
#'   zeros. Default = 0.8.
#' @param min_gene Minimal number of genes used in the initial phase. See
#'   details.
#' @param fix_num If true, always use \emph{min_gene} number of genes with the 
#'  highest proportion of non zeros in the initial phase. Default = F. See details.
#' @param cutoff The value below cutoff is treated as no expression. Default =
#'   0.1.
#' @param verbose Whether to show some intermediate results. Default = False.
#' @param num_mc The number of Gibbs steps for generating new imputed data after the
#'   parameters have been updated during Monte Carlo EM. Default = 3.
#' @param mcmc The number of Gibbs steps to sample imputed data after EM.
#'   Default = 50.
#' @param burnin The number of burnin steps before sample imputed data after EM.
#'   Default = 2.
#' @return \code{SIMPLE} returns a list of results in the following order.
#'   \enumerate{
#'     \item{loglik} {The log-likelihood of the full imputed gene expression at each iteration.}
#'     \item{pi} {The prior probabilites of cells belong to each cluster.}
#'     \item{mu} {Mean expression for each gene in each cluster}
#'     \item{sigma} {Variances of idiosyncratic noises for each gene in each cluster.}
#'     \item{beta} {Factor loadings.}
#'     \item{lambda} {Variances of factors for each cluster.}
#'     \item{z} {The posterior probability of each cell belonging to each cluster.}
#'     \item{Yimp0} {A matrix contains the expectation of gene
#'       expression specified by the model.}
#'     \item{pg} {A G by M0 matrix, dropout rate for each gene in each
#'     cluster estimated from initial clustering.}
#'     \item{initclus} {Output initial cluster results.}
#'     \item{impt} {A matrix contains the mean of each imputed
#'     entry by sampling multiple imputed values while the parameters are MLE. If mcmc <= 0, output
#'     the imputed expressoin matrix at last step of EM}
#'     \item{impt_var} {A matrix
#'     contains the variance of each imputed entry by sampling multiple imputed
#'     values while the parameters are MLE. NULL if mcmc <= 0.}
#'     \item{Ef} {If mcmc >0, output posterior means of factors
#'     given observed data (a n by K0 matrix). If mcmc <= 0, output conditional expectation of the factors for each cluster \eqn{E(f_i|z_i= m)} 
#'    at the last step of EM. A list with length M0, 
#'    each element in the list is a n by K0 matrix.}
#'     \item{Varf} {If mcmc >0, output posterior variances of
#'     factors given observed data (a n by K0 matrix). If mcmc <= 0, output conditional covariance matrix of factors for each cluster \eqn{Var(f_i|z_i = m)} at the last step of EM. 
#'      A list with length M0, each element in the list is a K0 by K0 matrix.}
#'     \item{consensus_cluster} {Score for the clustering stability of each cell by multiple imputations.
#'     NULL if mcmc <=0. }
#' }
#' @import doParallel
#' @importFrom foreach foreach
#' @seealso \code{\link{SIMPLE_B}}
#' @examples
#' library(foreach) 
#' library(doParallel) 
#' library(SIMPLEs) 
#'
#' # simulate number of clusters
#' M0 = 3 
#' # number of cells
#' n = 300
#' # simulation_bulk and getCluster is defined in the util.R under the util directory of the corresponding github repository.
#' source("utils/utils.R")
#' simu_data = simulation_bulk(n=300, S0 = 20, K = 6, MC=M0, block_size = 32, indepG = 1000 - 32*6, verbose=F, overlap=0)
#' Y2 = simu_data$Y2 
#' # number of factors
#' K0 = 6
#' # parallel
#' registerDoParallel(cores = 6)
#' # estimate the parameters and sample imputed values 
#' result <- SIMPLE(Y2, K0, M0, clus = NULL, K = 20, p_min = 0.5, max_lambda=T, min_gene = 200,cutoff=0.01)
#' # sample imputed values 
#' # evaluate cluster performance
#' celltype_true = simu_data$Z 
#' mclust::adjustedRandIndex(apply(result$z,1, which.max), celltype_true)
#' # or redo clustering based on imputed values (sometimes work better for real data)
#' getCluster(result$impt, celltype_true, Ks = 20, M0 = M0)[[1]]
#'
#' @author Zhirui Hu, \email{zhiruihu@g.harvard.edu}
#' @author Songpeng Zu, \email{songpengzu@g.harvard.edu}
#' @export

SIMPLE <- function(dat, K0, M0 = 1, iter = 10, est_lam = 1, impt_it = 5, penl = 1, 
    sigma0 = 100, pi_alpha = 1, beta = NULL, verbose = F, max_lambda = F, lambda = NULL, 
    sigma = NULL, mu = NULL, est_z = 1, clus = NULL, p_min = 0.8, cutoff = 0.1, K = 20, 
    min_gene = 300, num_mc = 3, fix_num = F, mcmc = 50, burnin = 2) {
    # EM algorithm initiation
    G <- nrow(dat)
    n <- ncol(dat)
    z <- NULL
    Y <- dat  # imputed matrix
    pg <- matrix(p_min, G, M0)
    
    pi <- rep(1/M0, M0)  # prob of z
    
    # random start
    if (is.null(lambda)) 
        lambda <- matrix(1/M0, M0, K0)  # sum to M0
    
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
        # fix number of hq genes for simulation, need to change back
        if (length(hq_ind) < min_gene) 
            hq_ind <- order(n1, decreasing = T)[1:min_gene]
    }
    
    # low dropout
    message(paste("inital impution for ", length(hq_ind), "high quality genes"))
    
    # init clustering
    if (is.null(clus)) {
        Y2_scale <- t(scale(t(dat[hq_ind, ])))
        s <- svd(Y2_scale)
        # for high dropout rate
        km0 <- kmeans(t(Y2_scale) %*% s$u[, 1:K], M0, iter.max = 80, nstart = 300)
        clus <- km0$cluster
    }
    
    z <- matrix(0, n, M0)
    for (m in 1:M0) z[clus == m, m] <- 1
    
    
    if (is.null(clus)) {
        res <- init_impute(dat[hq_ind, ], 1, rep(1, n), p_min, cutoff = cutoff, verbose = F)
        res[[2]] <- res[[2]] %*% t(rep(1, M0))
    } else {
        res <- init_impute(dat[hq_ind, ], M0, clus, p_min, cutoff = cutoff, verbose = F)
    }
    
    
    message("impute for hq genes")
    # iter, M0=1?
    impute_hq <- EM_impute(res[[1]], dat[hq_ind, ], res[[2]], M0, K0, cutoff, 20, 
        beta[hq_ind, ], sigma[hq_ind, , drop = F], lambda, pi, z, mu = NULL, celltype = clus, 
        penl, est_z, max_lambda, est_lam, impt_it, sigma0, pi_alpha, verbose = verbose, 
        num_mc = num_mc, lower = -Inf, upper = Inf)
    pg[hq_ind, ] <- res[[2]]
    beta[hq_ind, ] <- impute_hq$beta
    sigma[hq_ind, ] <- impute_hq$sigma
    # check this, as imputed mu may not be zero
    mu[hq_ind, ] <- impute_hq$mu
    Y[hq_ind, ] <- impute_hq$Y
    # gene_mean[hq_ind] = impute_hq$geneM
    z <- impute_hq$z
    
    
    nz <- colSums(impute_hq$z)
    Vm <- lapply(1:M0, function(m) impute_hq$Varf[[m]] * nz[m])
    
    # inital beta for other genes
    message("initial estimate beta for lq genes:")
    # which(n1 < p_min)
    lq_ind <- setdiff(1:G, hq_ind)
    # estimate beta and impute: only for positive part? (only impute for genes with
    # more than 10% nonzero) M step also estimate mu sapply(lq_ind, function(g){#
    res <- foreach(g = lq_ind, .combine = rbind) %dopar% {
        V <- 0
        for (m in 1:M0) V <- V + Vm[[m]]/sigma[g, m]
        W_temp <- c()
        Y_temp <- c()
        
        for (m in 1:M0) {
            Y_temp <- c(Y_temp, dat[g, ] * sqrt(z[, m])/sqrt(sigma[g, m]))
            Wb <- impute_hq$Ef[[m]] * sqrt(z[, m])/sqrt(sigma[g, m])
            Wmu <- matrix(0, n, M0)
            # n * M
            Wmu[, m] <- sqrt(z[, m])/sqrt(sigma[g, m])
            W_temp <- rbind(W_temp, cbind(Wmu, Wb))
        }
        
        ML <- cbind(matrix(0, K0, M0), chol(V))
        
        # (n+K) * (M + K)
        W_aug <- rbind(W_temp, ML)
        # G*(n+K)
        Y_aug <- c(Y_temp, rep(0, K0))
        
        
        # sigma^2
        penalty <- penl/(M0 * n + K0)/var(Y_aug)
        # K dimensional, n+K data
        fit1m <- glmnet(W_aug, Y_aug, family = "gaussian", alpha = 1, intercept = F, 
            standardize = F, nlambda = 1, lambda = penalty * K0/(M0 + K0), penalty.factor = c(rep(0, 
                M0), rep(1, K0)))
        
        coeff <- fit1m$beta[-1:-M0, 1]
        tempmu <- fit1m$beta[1:M0, 1]
        
        # sg = sum((Y_temp - fit1m$a0 - coeff %*% t(W_temp))^2) + sum(( coeff %*%V)*
        # coeff)* ns c(fit1m$a0, coeff, (sg+1)/(ns + 3))
        
        sg <- sapply(1:M0, function(m) {
            (sum((dat[g, ] - tempmu[m] - coeff %*% t(impute_hq$Ef[[m]]))^2 * z[, 
                m]) + sum((coeff %*% Vm[[m]]) * coeff))
        })
        c(fit1m$beta[, 1], rep((sum(sg) + 1)/(n + 3), M0))
    }
    
    mu[lq_ind, ] <- matrix(res[, 1:M0], ncol = M0)
    beta[lq_ind, ] <- matrix(res[, (M0 + 1):(K0 + M0)], ncol = K0)
    
    sigma[lq_ind, ] <- matrix(res[, -1:-(K0 + M0)], ncol = M0)
    sigma[sigma > 9] <- 9
    sigma[sigma < 1e-04] <- 1e-04
    
    
    # imputation set dropout rate as p_min
    if (M0 > 1) {
        # sample membership
        im <- apply(z, 1, function(x) which(rmultinom(1, 1, x) == 1))
    } else {
        im <- rep(1, n)
    }
    for (i in 1:n) {
        m <- im[i]
        ind <- which(dat[lq_ind, i] <= cutoff)
        ind <- lq_ind[ind]
        
        ms <- mu[ind, m] + beta[ind, , drop = F] %*% impute_hq$Ef[[m]][i, ]
        sds <- sqrt(sigma[ind, m])
        # need celltype
        p <- pg[ind, clus[i]]
        
        # compute x < 0 prob
        prob <- pnorm(cutoff, mean = ms, sd = sds)
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
        
        Y[ind, i] <- impt
    }
    
    message("EM for all genes")
    impute_result <- EM_impute(Y, dat, pg, M0, K0, cutoff, iter, beta, sigma, impute_hq$lambda, 
        impute_hq$pi, impute_hq$z, mu = NULL, celltype = clus, penl, est_z, max_lambda, 
        est_lam, impt_it = 1, sigma0, pi_alpha, verbose = verbose, num_mc = num_mc, 
        lower = -Inf, upper = Inf)
    
    impute <- matrix(0, n, G)
    for (m in 1:M0) {
        impute <- impute + t(impute_result$mu[, m] + impute_result$beta %*% t(impute_result$Ef[[m]])) * 
            impute_result$z[, m]
    }
    impute <- t(impute) * impute_result$geneSd + impute_result$geneM
    if (mcmc > 0) {
        message("multiple impution sampling")
        result2 = do_impute(dat, impute_result$Y, impute_result$beta, impute_result$lambda, 
            impute_result$sigma, impute_result$mu, impute_result$pi, impute_result$geneM, 
            impute_result$geneSd, clus, mcmc = mcmc, burnin = burnin, pg = pg, cutoff = cutoff, 
            verbose = verbose)

        return(list(loglik = impute_result$loglik, pi = impute_result$pi, mu = impute_result$mu, 
            sigma = impute_result$sigma, beta = impute_result$beta, lambda = impute_result$lambda, 
            z = impute_result$z, Yimp0 = impute, pg = pg, initclus = clus, impt = result2$impt, 
            impt_var = result2$impt_var, Ef = result2$EF, Varf = result2$varF, consensus_cluster = result2$consensus_cluster))
    }
    return(list(loglik = impute_result$loglik, pi = impute_result$pi, mu = impute_result$mu, 
        sigma = impute_result$sigma, beta = impute_result$beta, lambda = impute_result$lambda, 
        z = impute_result$z, Yimp0 = impute, pg = pg, initclus = clus, impt = impute_result$Y, 
        impt_var = NULL, Ef = impute_result$Ef, 
        Varf = impute_result$Varf, consensus_cluster = NULL))
}

