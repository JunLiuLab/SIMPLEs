#' Wrapper of SIMPLE or SIMPLE-B to select optimal K and M based on BIC


#' @param dat scRNASeq data matrix. Each row is a gene, each column is a cell.
#' @param bulk Bulk RNASeq data matrix. Should be log(1 + tpm/fpkm/rpkm). Each row is a gene which must be ordered the same as the scRNASeq data. Each column is a cell type which must be ordered as the cell type label in \emph{celltype}. Default: NULL
#' @param celltype A numeric vector for labels of cells in the scRNASeq. Each cell type has corresponding mean expression in the bulk RNASeq data. The labels must start from 1 to the number of types. If NULL, all cells are treated as a single cell type and the input bulk RNASeq should also have one column that is the mean expression over all the cell types. Default: NULL
#' @param b The scaling factor between scRNASeq and bulk RNASeq. If NULL, will do a weighted linear regression between mean of scRNASeq and bulk RNASeq for each cell type. Default = 1. Only relevant for SIMPLE-B. 
#' @param K0 a vector for number of latent gene modules to be selected.  Default is 10. Note that the minimum K0 is 2. 
#' @param M0 a vector for Number of clusters to be selected. Default is 1. If 1 is not in the sequence of M0, it will add 1 to the sequence; and the imputed matrix when M=1 is used to initialize imputation for other Ms. 
#' @param rel when increasing K, the algorithm will stop when the relative decrease of BIC is less than rel; but it will enumerate all given Ms. Default: 0.01.  
#' @param clus Initial clustering of scRNASeq data. If NULL, the function will
#'   use PCA and Kmeans to do clustering initially.
#' @param K The number of PCs used in the initial clustering. Default is 20.
#' @param iter The number of EM iterations using full data set. See details.
#' @param est_z The iteration starts to update Z.
#' @param impt_it The iteration starts to sample new imputed values in initial phase. See details.
#' @param max_lambda Whether to maximize over lambda. Default is True.
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
#'   zeros. Default = 0.4.
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


#' @return \code{selectKM} returns a matrix of BIC for all Ks and Ms tests, 
#'    mK, the best K; mM, the best M; result, 
#'    list for the result with smallest BIC in the following order.
#'   \enumerate{
#'     \item{loglik} {The log-likelihood of each MCMC sample of imputed gene expressionafter EM. NULL if mcmc <= 0}
#'     \item{loglik_tot} {The log-likelihood of the full imputed gene expression at each iteration and the prior of B matrix.}
#'     \item{BIC} {BIC which is -2 *loglik_tot + penalty on the number of parameters. Can be used to select paramters.}
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
#' @importFrom irlba irlba
#' @seealso \code{\link{SIMPLE}} \code{\link{SIMPLE_B}}

#' @examples
#' library(foreach)
#' library(doParallel)
#' library(SIMPLEs)
#'
#' # simulate number of clusters
#' M0 <- 3
#' # number of cells
#' n <- 300
#' # simulation_bulk and getCluster is defined in the util.R under the util directory of the corresponding github repository.
#' source("utils/utils.R")
#' simu_data <- simulation_bulk(n = 300, S0 = 20, K = 6, MC = M0, block_size = 32, indepG = 1000 - 32 * 6, verbose = F, overlap = 0)
#' Y2 <- simu_data$Y2
#' # number of factors
#' K <- c(6,10)
#' M <- c(1, 3)
#' # parallel
#' registerDoParallel(cores = 6)
#' # estimate the parameters and sample imputed values
#' results <- selectKM(Y2, K, M, clus = NULL, K = 20, p_min = 0.5, max_lambda = T, min_gene = 200, cutoff = 0.01)
#' print(sprintf("best M and K: %d, %d", results$mM, results$mK))
#' result = results$result
#' # evaluate cluster performance
#' celltype_true <- simu_data$Z
#' mclust::adjustedRandIndex(apply(result$z, 1, which.max), celltype_true)
#' # or redo clustering based on imputed values (sometimes work better for real data)
#' getCluster(result$impt, celltype_true, Ks = 20, M0 = M0)[[1]]
#' @author Zhirui Hu, \email{zhiruihu@g.harvard.edu}
#' @author Songpeng Zu, \email{songpengzu@g.harvard.edu}

#' @export
selectKM <- function(dat, bulk = NULL, celltype = NULL, b = 1,  K0 = 10, M0 = 1, iter = 10, est_lam = 1, impt_it = 5, penl = 1,
                   sigma0 = 100, pi_alpha = 1, beta = NULL, verbose = F, max_lambda = T, lambda = NULL,
                   sigma = NULL, mu = NULL, est_z = 1, clus = NULL, p_min = 0.4, cutoff = 0.1, K = 20,
                   min_gene = 300, num_mc = 3, fix_num = F, mcmc = 50, burnin = 2, rel = 0.01) {

	G <- nrow(dat)
	n <- ncol(dat)
	mK = mM = NULL
  mBIC = Inf
  best = list() # for each M
  M0 = sort(M0)
  if(!(1 %in% M0)) M0 = c(1, M0)
  K0 = sort(K0)

  #record all BIC
  BICs = data.frame(matrix(NA, nrow = length(M0), ncol = length(K0)))
  rownames(BICs) = as.character(M0); colnames(BICs) = as.character(K0)

  init_imp = NULL
  for(M1 in M0)
  {
     prevBIC = Inf
     print(sprintf("number of clusters: %d", M1))
     for(K1 in K0)
     {
       print(sprintf("number of factors: %d", K1))
       result = SIMPLE(dat, K1, M1, iter = iter, est_lam = est_lam, impt_it = impt_it, penl = penl, init_imp = init_imp, 
           sigma0 = sigma0, pi_alpha = pi_alpha, beta = beta, verbose = verbose, max_lambda = max_lambda, lambda = lambda,
           sigma = sigma, mu = mu, est_z = est_z, clus = clus, p_min = p_min, cutoff = cutoff, K = K,
           min_gene = min_gene, num_mc = num_mc, fix_num = fix_num, mcmc = mcmc, burnin = burnin) 
       if(is.null(p_min)) p_min = result$p_min
       print(sprintf("BIC: %.0f", result$BIC0))
       if(n < 1000)
       {
          if(mBIC > result$BIC0)
          {
             mBIC = result$BIC0
             mK = K1
             mM = M1 
             best[[M1]] = result
          }

          BICs[as.character(M1), as.character(K1)] = result$BIC0
          if((result$BIC0 - prevBIC) > -abs(result$BIC0)* rel) break;
          prevBIC = result$BIC0
        }else{
         if(mBIC > result$BIC)
         {
             mBIC = result$BIC
             mK = K1
             mM = M1 
             best[[M1]] = result
          }
          BICs[as.character(M1), as.character(K1)] = result$BIC
          if((result$BIC - prevBIC) > -abs(result$BIC)* rel) break;
          prevBIC = result$BIC
        }
		}

    if(M1 == 1) init_imp = best[[1]]$impt
	}

  return(list(BICs = BICs, mK = mK, mM = mM, result = best[[mM]]))
}


