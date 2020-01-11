
# fit-truncated normal given sigma solve mu
LL <- function(par, xbar, Sx, delta = 0) {
    sigma <- par[1]
    
    mu <- xbar + (Sx - sigma^2)/(xbar - delta)
    
    return(log(sigma) + (Sx + (xbar - mu)^2)/2/sigma^2 + pnorm((mu - delta)/sigma, 
        log.p = T))
}

# LL1 <- function(par, xbar, Sx, r01, p, delta = 0) # { mu = par[1] sigma =
# par[2] phi = pnorm( (delta-mu)/sigma) log(sigma) + (Sx + (xbar -
# mu)^2)/2/sigma^2 - r01*log(1 - p + p*phi) }

# fit zero-inflated censored normal with p fixed (at either p_min or p_max),
# given sigma solve mu
LLwZ <- function(par, xbar, Sx, r01, p, delta = 0) {
    sigma <- par[1]
    
    mu <- xbar + (Sx - sigma^2)/(xbar - delta)
    
    phi <- pnorm((delta - mu)/sigma)
    log(sigma) + (Sx + (xbar - mu)^2)/2/sigma^2 - r01 * log(1 - p + p * phi)
}

#' Fit a zero-inflated censored normal distribution.
#'
#' Fit each gene expression by a zero-inflated censored normal distribution by
#' direct maximize the log-likelihood and return MLE of mean, variance and
#' probability of mixture component.
#'
#' @details
#' It will return NA if number of nonzero is less than 4 indicating that the gene has no expression
#' or if the optimization is not convergence or error
#' which will be filled with sample mean and std, p_max afterwards.
#'
#' @param dat Data matrix.
#' @param cutoff The value below cutoff is treated as no expression. Default =
#'   0.1.
#' @param iter Max number of iteration for optimization. Default = 100.
#' @param s_upper upper bound of standard deviation. Default = 4.
#' @param s_lower lower bound of standard deviation. Default = 0.1.
#' @param p_min lower bound of the probability of the normal component. Default
#'   = 0.6.
#' @param p_max upper bound of the probability of the normal component. Can be either a scaler or a vector. If it is a vector, each element is the upper bound for a gene. Default = 1.
#' @return \code{ztruncnorm} returns a list of results in the following order.
#' \enumerate{
#'   \item{param}{A matrix with mean and standard deviation for each gene. }
#'   \item{pg}{A vector for the probablity of the normal component.}
#' }
#'
#' @author Zhirui Hu, \email{zhiruihu@g.harvard.edu}
#' @author Songpeng Zu, \email{songpengzu@g.harvard.edu}
# p_max is the dropout rate for integrating bulk
ztruncnorm <- function(dat, cutoff = 0.1, iter = 200, s_upper = 4, s_lower = 0.1, 
    p_min = 0.6, p_max = 1) {
    n <- ncol(dat)
    dat_pos <- dat
    dat_pos[dat_pos <= cutoff] <- NA
    mu <- rowMeans(dat_pos, na.rm = T)
    var <- rowVars(dat_pos, na.rm = T)
    n1 <- rowSums(!is.na(dat_pos))
    var <- var * (n1 - 1)/n1
    
    if (length(p_max) == 1) 
        p_max <- rep(p_max, nrow(dat))
    
    # first fit truncated normal distribution using only positive entries write the
    # log-likelihood in terms of only std and optimize by Brent
    param0 <- foreach(g = 1:nrow(dat_pos), .combine = rbind) %dopar% {
        # if(is.na(mu[g]) | is.na(var[g])) return(c(NA,NA))
        if (n1[g] < 4) 
            {
                return(c(NA, NA))
            }  # return(c(0,cutoff/2))
        # if(mu[g] - cutoff > 3.5*sqrt(var[g])) return(c(mu[g], sqrt(var[g])))
        tryCatch({
            res <- optim(sqrt(var[g]), LL, xbar = mu[g], Sx = var[g], delta = cutoff, 
                method = "Brent", upper = s_upper, lower = s_lower, control = list(maxit = iter))
            if (res$convergence != 0) {
                return(c(NA, NA))
            } else {
                c(mu[g] + (var[g] - res$par[1]^2)/(mu[g] - cutoff), res$par[1])
            }
        }, error = function(err) {
            c(NA, NA)
        })
    }
    
    # refit if estimated p_g is too large or too small
    n0 <- n - n1
    r01 <- n0/n1
    p1 <- n1/n
    p11 <- pnorm((param0[, 1] - cutoff)/param0[, 2])  # is na when too few nonzero
    pg <- p1/p11  ## 1 - (p0 - p01)
    
    ind <- which(pg > p_max | is.na(pg))  # p11 == na or p1 = p11 = 0
    pg[ind] <- p_max[ind]
    ind2 <- which(pg < p_min)
    pg[ind2] <- p_min
    param1 <- foreach(g = c(ind, ind2), .combine = rbind) %dopar% {
        if (is.na(var[g])) 
            {
                return(c(-1, 0.5))
            }  # cutoff/2
        tryCatch({
            res <- optim(par = sqrt(var[g]), fn = LLwZ, xbar = mu[g], Sx = var[g], 
                r01 = r01[g], p = pg[g], delta = cutoff, method = "Brent", upper = s_upper, 
                lower = s_lower, control = list(maxit = iter))
            # res <- optim(par=c(mu[g], sqrt(var[g])), fn=LL1, xbar = mu[g], Sx = var[g], r01
            # = r01[g], p = pg[g], delta = cutoff, control = list(maxit = iter))
            if (res$convergence != 0) {
                return(c(NA, NA))
            } else {
                c(mu[g] + (var[g] - res$par[1]^2)/(mu[g] - cutoff), res$par[1])
                # res$par
            }
        }, error = function(err) {
            c(NA, NA)
        })
    }
    
    param0[c(ind, ind2), ] <- param1
    
    return(list(param = param0, pg = pg))
}
