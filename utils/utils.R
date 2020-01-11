# need ggplot
simulation_bulk <- function(n = 300, S0 = 10, K = 3, MC = 2, block_size = 50, overlap = 15,
    indepG = 30, dropr = 0.3, verbose = F) {
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
    
    # Sigma = matrix(0, G, G) diag(Sigma) <- rnorm(G, 0.6, 0.1)^2
    Sigma = rgamma(G, 2, 1/0.3)  #rnorm(G, 0.6, 0.1)
    
    # Mu = rtnorm(G, lower = 0, mean = 1.8, sd = 0.5)%*% t(rep(1,MC))
    Mu = exp(rnorm(G, mean = 0.5, sd = 0.5)) %*% t(rep(1, MC))
    
    act_ind = sample(1:G, S0 * MC)
    label0 = matrix(0, G, MC)
    # specify mean for each cluster
    for (m in 1:MC) {
        ss = ((m - 1) * S0 + 1):(m * S0)
        Mu[act_ind[ss], m] = Mu[act_ind[ss], m] * sample(c(0.2, 0.5, 1.2, 1.5, 2), 
            S0, replace = T)  #c(0.3, 0.5, 2 , 3, 4)
        label0[act_ind[ss], m] = 1
    }
    
    
    if (K > 0) {
        # Factor Loading
        Gamma <- matrix(0, G, K)
        Gamma[1:block_size, 1] <- 1
        if (K > 1) {
            for (i in 2:K) {
                Gamma[((block_size - overlap) * (i - 1) + 1):(block_size + (block_size - 
                  overlap) * (i - 1)), i] <- 1
            }
        }
        
        B = Gamma/4  #sqrt(block_size) # eigenvalue of BB^T = block_size * B^2
        
        Lambda = list()  # specify variance for each cluster
        for (m in 1:MC) {
            Lambda[[m]] = diag(1, K, K)
        }
        
        W = sapply(Z, function(z) mixtools::rmvnorm(1, rep(0, K), Lambda[[z]]))  #matrix(rnorm(K*n), K, n) # K * n
        E = matrix(rnorm(G * n), nrow = G) * Sigma  #t(mvrnorm(n,rep(0,G),Sigma))
        Y = Mu[, Z] + B %*% W + E
        
    } else {
        Y = matrix(rnorm(G * n, Mu[, Z], Sigma), nrow = G)
    }
    
    
    # add dropout
    Y2 = Y
    dZ = matrix(rbinom(G * n, 1, exp(-dropr * rowMeans(Mu)^2)), nrow = G)  # 0.1
    Y2[dZ == 1] = 0
    Y2[Y2 < 0] = 0
    
    ind = which(rowSums(Y2 != 0) > 4)
    Y2 = Y2[ind, ]
    Y = Y[ind, ]
    label0 = label0[ind, ]
    B = B[ind, ]
    Mu = Mu[ind, ]
    Sigma = Sigma[ind]
    
    # show Gamma and PCA of Y
    if (verbose) {
        if (K > 0) 
            plot(ggplot(melt(Gamma), aes(Var2, Var1)) + geom_tile(aes(fill = value)) + 
                scale_fill_gradient(low = "white", high = "red", guide = "colourbar") + 
                theme_bw() + xlab("") + ylab(""))
        svds = svd(t(scale(t(Y2))), nu = 20, nv = 20)
        par(mfrow = c(2, 2))
        plot(svds$d[1:20])
        plot(svds$v[, 1:2], col = Z, pch = 16)
        hist(rowMeans(Y2 > 0))
    }
    
    return(list(Y2 = Y2, Y = Y, B = B, W = W, Mu = Mu, Lambda = Lambda, Sigma = Sigma, 
        Z = Z, bulk = rowMeans(Mu[, Z]), S_label = label0))  #
}

getCluster <- function(Y2, celltype, Ks = NULL, M0 = NULL) {
    if (is.null(M0)) 
        M0 = length(unique(celltype))
    if (is.null(Ks)) 
        Ks = M0 - 1
    Y2_scale = t(scale(t(Y2)))
    s = svd(Y2_scale)
    # hist(rowSds(Y2_scale))
    result <- matrix(0, length(Ks), 2)
    i = 1
    for (K in Ks) {
        km <- kmeans(t(Y2_scale) %*% s$u[, 1:K], M0, iter.max = 80, nstart = 300)
        km0 <- kmeans(s$v[, 1:K], M0, iter.max = 80, nstart = 300)
        result[i, ] = c(mclust::adjustedRandIndex(km$cluster, celltype), mclust::adjustedRandIndex(km0$cluster, 
            celltype))
        i = i + 1
    }
    
    
    return(list(result, xtabs(~km0$cluster + celltype)))
}

# tsne plot n*G, K
tsneplot <- function(x, celltype, pca = T, scale = T, center = T, K = 20, file = NULL) {
    
    result <- Rtsne::Rtsne(x, pca_scale = scale, pca_center = center, initial_dims = K, 
        pca = pca)
    dat = data.frame(cbind(result$Y[, 1:2]), celltype)
    colnames(dat) = c("TSNE1", "TSNE2", "Type")
    p <- ggplot(dat, aes(x = TSNE1, y = TSNE2, color = Type)) + geom_point() + theme_bw()
    
    if (!is.null(file)) {
        pdf(file)
        plot(p)
        dev.off()
    } else {
        plot(p)
    }
    
}

# modified t-test
myttest <- function(x, m, celltype) if (sd(x) < 1e-06) {
    1
} else {
    t.test(x ~ celltype == m)$p.value
}

# get mse of correla
mycor <- function(X, m, celltype, cors) {
    org_cors = cor(t(X[, celltype == m]))
    org_cors[is.na(org_cors)] = 0
    return(c(mean(org_cors[cors == 0]^2), mean((org_cors[(cors != 0) & (cors != 1)] - 
        cors[(cors != 0) & (cors != 1)])^2)))
}

mycor2 <- function(X, m, celltype, cors) {
    org_cors = cor(t(X[, celltype == m]))
    org_cors[is.na(org_cors)] = 0
    return(mean((org_cors[(cors != 1)] - cors[(cors != 1)])^2, na.rm = T))
}
