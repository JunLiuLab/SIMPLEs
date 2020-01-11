# simulation according to scIMpute
library(scImpute)
library(Rmagic)
library(foreach)
library(doParallel)
library(SIMPLE)

source("Simulation/utils.R")

## preprocessing ##
dat <- read.csv("chu_data/GSE75748_sc_cell_type_ec.csv", row.names = 1)
bulk <- read.csv("chu_data/GSE75748_bulk_cell_type_ec.csv",row.names = 1)

# remove ERCC and MT
er_id0 = grep("^ERCC", rownames(bulk))
er_id = grep("^ERCC", rownames(dat))

mt_id0 =  grep("^MT-", rownames(bulk))
mt_id =  grep("^MT-", rownames(dat))


dat = dat[-c(er_id, mt_id),]
bulk = bulk[-c(er_id0, mt_id0),]

ind = which(rowMeans(dat>0) > 0.1)
dat = dat[ind, ]
write.csv(dat, "dat_filter.csv")

# get differential expressed genes in bulk
library(DESeq2)
condition = factor(sapply(colnames(bulk), function(x) strsplit(x,"_")[[1]][1]))
#condition = plyr::mapvalues(condition, from="H9", to="H1")
dds <- DESeqDataSetFromMatrix(countData = round(bulk),
                              colData = DataFrame(condition),
                              design= ~ condition)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="condition_EC_vs_DEC")

cl = 7
true_DE = matrix(0, nrow(dat), cl) #bulk
rownames(true_DE) = rownames(dat) #bulk

for(i in 1:cl)
{
  ct = c(0,rep(-1/(cl-1),(cl-1)))
  if(i>1) ct[i] = 1

  res <- results(dds, contrast=ct)

  # select uprelated genes
  res <- res[rownames(dat), ]
  true_DE[,i] = res$padj * (2* (res$log2FoldChange >0) - 1)
  true_DE[abs(res$log2FoldChange) < 0.5,i] = 1

}
colnames(true_DE) = levels(condition)
write.csv(true_DE, file="chu_data/bulk_diff.csv")


dat_norm <- t(log(t(dat)/colSums(dat)*1e6+1))
bulk_norm = t(log(t(bulk)/colSums(bulk)*1e6+1))


ind = which(rowMeans(dat_norm>0) > 0.1)
print(length(ind))


dat_norm = dat_norm[ind, ]


celltype= sapply(colnames(dat_norm), function(x) strsplit(x,"_")[[1]][1])
##batch =sapply(colnames(dat_norm), function(x) str_extract(x,"(Exp|Batch)\\d+")[[1]][1])

celltype0 = sapply(colnames(bulk_norm), function(x) strsplit(x,"_")[[1]][1])

bulk_norm = bulk_norm[rownames(dat_norm),]
bulk_norm_mean = sapply(levels(factor(celltype)), function(x) rowMeans(bulk_norm[, celltype0==x]))

save(bulk_norm, dat_norm, celltype, celltype0, bulk_norm_mean, file="chu_data/data.rdat")

# merge with bulk
dat_join = merge(dat_norm, bulk_norm, by=0)
rownames(dat_join) = dat_join[,1]
dat_join = dat_join[,-1]

for(x in levels(factor(celltype)))
{
  temp_dat = data.frame("scp" = rowMeans(dat_pos[,celltype==x],na.rm = T), "sc" = rowMeans(dat_norm[,celltype==x]), "bulk" = rowMeans(bulk_norm[rownames(dat_norm),celltype0==x]), "impute" = rowMeans(res[[1]][,celltype==x]))

  ggplot(temp_dat, aes(x=bulk, y = impute)) + geom_point() + stat_smooth(aes(x=bulk, y = impute))

  ggplot(temp_dat, aes(x=sc, y = impute)) + geom_point() + stat_smooth(aes(x=sc, y = impute))

  print(cor(temp_dat$impute,temp_dat$bulk))
  print(cor(temp_dat$sc,temp_dat$bulk))
  print(cor(temp_dat$scp,temp_dat$bulk, use="pairwise.complete.obs"))

}

par(mfrow=c(3,3), mar=c(3,3,2,1))
for(i in 1:length(er_id))
{

  #boxplot(unlist(dat_norm[er_id[i],])~celltype, main=rownames(dat_norm)[er_id[i]])
  boxplot(unlist(bulk_norm[er_id0[i],])~celltype0, main=rownames(bulk_norm)[er_id0[i]])
}
boxplot(colSums(dat)~celltype)
boxplot(colSums(dat_norm>0)~celltype)




## tsne using original data
#Y2_scale = t(scale(t(dat_norm[, ])))
Y2_scale = t(scale(t(Y2)))
s = svd(Y2_scale)
s = svd(result$Y)
s = svd(sci_result[[1]])
ts = Rtsne::Rtsne(s$v[,1:20])

clus0 = apply(result$z, 1, which.max)
pdf("factor_model_result/chu_scimpute_K6_tsne.pdf") #
ggplot(data.frame(ts$Y[,1:2], 'celltype' = factor(celltype_true)), aes(X1,X2,color=celltype)) + geom_point() + theme_bw() + xlab("tSNE1") + ylab("tSNE2")
ggplot(data.frame(ts$Y[,1:2], 'celltype' = factor(clus0)), aes(X1,X2,color=celltype)) + geom_point() + theme_bw() + xlab("PC1") + ylab("PC2") #s$v[,1:2]
dev.off()

getCluster(Y2[sel,], celltype2, Ks=c(10,15,20,25))


tsneplot(t(Y2), celltype2, scale=T, center=T, K=20)
tsneplot(result$Ef[[1]], celltype2, scale=F, center=F)
tsneplot(sci_result[[1]], celltype2, scale=T, center=T, K=20)




#### step1: impute and estimate B ####
run = 5
cluster_result <- matrix(0, run, 8)
mse <- matrix(0, run, 4)
auc <- matrix(0, run, 4)


label = matrix(0, G, M0)

for(m in 1:M0)
{
  ss = ((m-1) * S0 + 1) : (m * S0)
  label[ss ,m] = 1
}

## clustering performance
getCluster(res[[1]][step1[[2]],], celltype, Ks=c(2:10,20))


temp_data = data.frame("bulk" = bulk_norm_mean[,1],"sc" = rowMeans(res[[1]][,celltype=="DEC"]), "d" = res[[2]][,1])


n = ncol(dat_norm)
prop = xtabs(~celltype)/n
bulk = bulk_norm_mean %*% prop

ggplot(temp_data, aes(x=bulk, y = sc, color = d)) + geom_point() + scale_color_gradient2()

load("../chu_data/data.rdat")
K0 = 10
M0 = 7
pm = 0.7

# combine H1 and H9 cell types
celltype0 = celltype
celltype = plyr::mapvalues(celltype, from="H9", to="H1")
bulk_norm_mean0 = bulk_norm_mean
bulk_norm_mean[,3] = (bulk_norm_mean[,3] + bulk_norm_mean[,4])/2
bulk_norm_mean = bulk_norm_mean[,-4]
# t-test using orignial data set as true differential genes
i=1
true_DE = matrix(0, nrow(dat_norm), 6)

for(m in levels(factor(celltype)))
{
  ind0 = which(celltype==m)
  pvs = apply(dat_norm, 1, function(x) {
    if(mean(x[ind0]) < mean(x[-ind0]) + 0.2) return(1)
    wilcox.test(x[ind0], x[-ind0], alternative="g")$p.value #t
    })
  fdr = p.adjust(pvs)
  true_DE[fdr<0.05 ,i] = 1
  i = i+1

}


for( r in 1:run)
{
  # add more dropouts to original data and subsample cells
  cell_ind = sample(1:ncol(dat_norm), ncol(dat_norm)/2)
  Y2 = dat_norm[, cell_ind]
  celltype_true = celltype[cell_ind]
  celltype20 = celltype0[cell_ind]

  # # add dropout
  # G = nrow(Y2)
  # n = ncol(Y2)
  # Z = matrix(rbinom(G*n,1,exp(-0.5*Y2)), nrow=G) #0.3 * exp,0.5
  # #Z = matrix(rbinom(G*n,1,0.2), nrow=G)
  # Y2[Z==1] = 0

  #ind = which(rowSums(Y2 != 0) > 50)
  # random sample 2000 genes
  ind = sample(1:nrow(dat_norm), 3000)
  Y2 = Y2[ind, ]
  Y = dat_norm[ind, cell_ind]
  #Z2 = Z[ind, ]
  bulk_mean = bulk_norm_mean[ind,]

  n = ncol(Y2)
  prop = xtabs(~celltype_true)/n
  bulk = bulk_mean %*% prop

  # # get bulk
  # prop = xtabs(~celltype2)/n
  # bulk = bulk_norm_mean %*% prop
  # bulk = bulk[ind,,drop=F ]
  # DE = true_DE[ind,]

  # registerDoParallel(6)
  # result <- ztruncnorm0(Y2, cutoff = 0)
  #
  # step1 <- IF(Y2, data.frame(result), norm =T, colcenter = F, opt = 2, verbose=T, cutoff = 0)

  sel = which(rowSds(Y2)>1.5)
  Y2 = Y2[sel, ]
  Y = Y[sel, ]
  Z2 = Z2[sel, ]
  #bulk = bulk[sel,,drop=F ]
  bulk_mean = bulk_mean[sel,]
  #DE = DE[sel,]

  #### method 4 ####
  registerDoParallel(cores = 7)
  result <- scimpclu(Y2, K0, M0, clus = NULL, K = 20, iter= 10, est_z = 1, impt_it = 1, max_lambda=T, est_lam = 1, penl =1, sigma0 = 100, p_min = pm, cutoff=0.1,verbose=T, min_gene = 300)
  resultY = result$Y

  resultB <- scimpclu_bulk(Y2, K0, bulk, M0, celltype=rep(1, n), clus = NULL, K = 20, iter= 10, est_z = 1, impt_it = 1, max_lambda=T, est_lam = 2, penl =1, sigma0 = 100, p_min = pm, min_gene = 300,cutoff=0.1,verbose=T)
  resultYB = resultB$Y


  # result <- doEM_impute(Y2, 6, 20, celltype = celltype2, iter=20, est_lam = 3, verbose=T, max_lambda=F, penl =2, clus = NULL, est_z = 1,sigma0 = 100, p_min =0.7, cutoff = 0.1, K = 20) #6 penl=5
  #
  #
  # tsneplot(result$Ef[[1]])
  #
  # result <- doEM_impute_bulk(Y2, 7, K0, bulk_mean, as.numeric(factor(celltype2)), iter=0, est_lam = 3, verbose=T, max_lambda=F, penl =1, clus = NULL, est_z = 1,sigma0 = 100, p_min =0.7) #dat_norm[step1[[2]],]
  #
  # ## step2: impute ##
  # # result2 = do_impute(Y2, result$Y, result$posM, result$posSd, 6, K0, result$beta, result$lambda, result$sigma, result$mu, result$pi, pg = result$pg, mcmc=90, burnin = 10, verbose=T, cutoff = 0) #M0, K0,
  #
  # result2 = do_impute(Y2, result$Y, result$geneM, result$geneSd, 1, 20, result$beta, result$lambda, result$sigma, result$mu, result$pi, pg = result$pg, mcmc=90, burnin = 10, verbose=T, cutoff = 0.1) #result$pg, 0.8
  #
  # ## with bulk
  # result2 = do_impute(Y2, result$Y, result$geneM, result$geneSd, 1, K0, result$beta, result$lambda, result$sigma, result$mu, result$pi, mcmc=90, burnin = 10, verbose=F, cutoff = 0.5, pg=result$pg, clus = as.numeric(factor(celltype2))) #M0



  #### method 3 ####
  sci_result <- scimpute("", infile = "csv", outfile = "csv", type = "count",
                         "../scimpute/", labeled = FALSE, drop_thre = pm, Kcluster = 7,
                         labels = NULL, genelen = NULL, ncores =6, count_lnorm = Y2+log10(1.001))#7


  getCluster(sci_result[[1]], celltype_true, Ks = c(10,20,30), M0 = 7)[[1]]


  #sci_result2 <- scimpute("", infile = "csv", outfile = "csv", type = "count",
                         # "scimpute/", labeled = FALSE, drop_thre = 0.8, Kcluster = 1,
                         # labels = NULL, genelen = NULL, ncores =6, count_lnorm = Y2+log10(1.001))
  #save(dat_norm, step1, sci_result, file="chu_data/scimpute1.rdat")

  save(Y2, Z2, sci_result, file="chu_data/scimpute2.rdat")
  save(Y2, sci_result, file="chu_data/scimpute_part_cells_K6.rdat")
  save(Y2, Y, Z2, cell_ind, result, result2, file="chu_data/fa_impute_part_cells_K6_bulk.rdat")

  #sci_result <- read.csv("chu_data/scimpute_count.csv", row.names = 1)

  #### method 2 ####
  MAGIC <- magic(t(Y2), genes="all_genes", t="auto")

  ## cluster (project onto 50 dim space)
  getCluster(Y2, celltype2, Ks=c(15,20,25,30))
  getCluster(sci_result[[1]], celltype2, Ks=c(5, 10,15, 20))
  getCluster(result2$impt, celltype2, Ks=c(5, 10,15, 20))

  cluster_result[r,7:8] <- colMaxs(getCluster(result2$impt, celltype2, Ks=c(20,25,28,30,40)))
  cluster_result[r,5:6] <- colMaxs(getCluster(sci_result[[1]], celltype2, Ks=c(5,10,15,20,30)))
  cluster_result[r,3:4] <- colMaxs(getCluster(t(MAGIC$result), celltype2, Ks=c(5,10,15,20,25)))
  cluster_result[r,1:2] <- colMaxs(getCluster(Y2, celltype2, Ks=c(10,15,20,25)))

  # correlation within clusters
  cor0 = sapply(unique(celltype2), function(x) {
    a = cor(Y[, celltype2 == x])
    mean(a[upper.tri(a)])
  }
  )

  cor1 = sapply(unique(celltype2), function(x) {
    a = cor(Y2[, celltype2 == x])
    mean(a[upper.tri(a)])
  }
  )

  cor2 = sapply(unique(celltype2), function(x) {
    a = cor(t(MAGIC$result)[, celltype2 == x])
    mean(a[upper.tri(a)])
  }
  )

  cor3 = sapply(unique(celltype2), function(x) {
    a = cor(sci_result[[1]][, celltype2 == x])
    mean(a[upper.tri(a)])
  }
  )

  cor4 = sapply(unique(celltype2), function(x) {
      a = cor(result2$impt[, celltype2 == x])
      mean(a[upper.tri(a)])
    }
    )

  # correlation between clusters
  cl = levels(factor(celltype2))
  cor_bet = rep(0, 5)
  for(m in 1:6)
  {
    ind_m = which(celltype2 == cl[m])
    for(m2 in (m+1):7)
    {
        ind_m2 = which(celltype2 == cl[m2])
        cor_bet[1] = cor_bet[1] + mean(cor(Y[, ind_m], Y[, ind_m2]))
        cor_bet[2] = cor_bet[2] + mean(cor(Y2[, ind_m], Y2[, ind_m2]))
        cor_bet[3] = cor_bet[3] + mean(cor(t(MAGIC$result)[, ind_m], t(MAGIC$result)[, ind_m2]))
        cor_bet[4] = cor_bet[4] + mean(cor(sci_result[[1]][, ind_m], sci_result[[1]][, ind_m2]))
        cor_bet[5] = cor_bet[5] + mean(cor(result2$impt[, ind_m], result2$impt[, ind_m2]))
    }
  }


 # gene-gene correlation


  # correlation between sc and bulk
  corB = matrix(0,6,5)
  i = 1
  for( x in levels(factor(celltype2)))
  {
    print(x)
    corB[i,4] = cor(bulk_mean[,i], rowMeans(sci_result[[1]][,celltype2 == x]))
    corB[i,3] = cor(bulk_mean[,i], colMeans(MAGIC$result[celltype2 == x,]))
    corB[i,5] = cor(bulk_mean[,i], rowMeans(result2$impt[, celltype2 == x]))
    corB[i,1] = cor(bulk_mean[,i], rowMeans(Y[, celltype2 == x]))
    corB[i,2] = cor(bulk_mean[,i], rowMeans(Y2[, celltype2 == x]))
    i = i+1
  }

  #mse
  mean(Y[Z2==1])
  mean((t(MAGIC$result)[Z2==1] - Y[Z2==1])^2)
  mean(abs(sci_result[[1]][Z2==1 & Y>0] - Y[Z2==1 & Y>0]))
  imt = result2$impt
  imt[imt <0] = 0
  mean(abs(imt[Z2==1 & Y>0] - Y[Z2==1 & Y>0]))

  mean((imt[Z2==1] - Y[Z2==1])^2)


  zz = which(Z2==1 & Y==0)
  nz = which(Z2==1 & Y!=0)

  # auc
  auc = matrix(0,5,5)
  i = 1
  for(m in levels(factor(celltype2)))
  {
    ind0 = which(celltype2==m)

    fdr0 = p.adjust(apply(Y, 1, function(x) {
      if(mean(x[ind0]) < mean(x[-ind0] + 0.1) ) return(1)#+ 0.2
      #t.test(x[ind0], x[-ind0], alternative="g")$p.value
      wilcox.test(x[ind0], x[-ind0], alternative="g")$p.value
    }))

    fdr1 = p.adjust(apply(Y2, 1, function(x) {
      if(mean(x[ind0]) < mean(x[-ind0]+ 0.1)) return(1)
      wilcox.test(x[ind0], x[-ind0], alternative="g")$p.value
    }))

    fdr4 = p.adjust(apply(result2$impt, 1, function(x) {
      if(mean(x[ind0]) < mean(x[-ind0])+ 0.1 ) return(1)
      wilcox.test(x[ind0], x[-ind0], alternative="g")$p.value
    }))

    # fdr2 = p.adjust(apply(result2$impt, 1, function(x) {
    #   if(mean(x[ind0]) < mean(x[-ind0]) + 0.1) return(1)
    #   t.test(x[ind0], x[-ind0], alternative="g")$p.value
    # }))

    fdr2 = p.adjust(apply(t(MAGIC$result), 1, function(x) {
      if(mean(x[ind0]) < mean(x[-ind0])+ 0.1 ) return(1)
      wilcox.test(x[ind0], x[-ind0], alternative="g")$p.value
    }))

    fdr3 = p.adjust(apply(sci_result[[1]], 1, function(x) {
      if(mean(x[ind0]) < mean(x[-ind0])+ 0.1) return(1)
      wilcox.test(x[ind0], x[-ind0], alternative="g")$p.value
    }))

    auc[r, 1] = auc[r, 1] + caTools::colAUC(fdr0, DE[,i])
    auc[r, 2] = auc[r, 2] + caTools::colAUC(fdr1, DE[,i])
    auc[r, 3] = auc[r, 3] + caTools::colAUC(fdr2, DE[,i])
    auc[r, 4] = auc[r, 4] + caTools::colAUC(fdr3, DE[,i])
    auc[r, 5] = auc[r, 5] + caTools::colAUC(fdr4, DE[,i])

    i = i+1

  }


  # only compare DEC and EC
  ind0 = which(celltype == "DEC")
  ind1 = which(celltype == "EC")

  fdr_t = p.adjust(apply(dat_norm[rownames(Y),], 1, function(x) {
    if(mean(x[ind0]) < mean(x[ind1])+ 0.2) return(1)
    wilcox.test(x[ind0], x[ind1], alternative="g")$p.value
  }))

  DE_t = fdr_t<0.05 #[ind[sel]]

  # select uprelated genes
  res <- results(dds, name="condition_EC_vs_DEC")
  res <- res[rownames(Y2), ]
  #res <- res[sel,]
  DE_t = rep(0, nrow(res))
  DE_t[which(res$log2FoldChange < -1 & res$padj <1e-3)] = 1


  ind0 = which(celltype2 == "DEC")
  ind1 = which(celltype2 == "EC")
  pval1 = apply(result2$impt, 1, function(x) {
    if(mean(x[ind0]) < mean(x[ind1])+0.05 ) return(1)
    wilcox.test(x[ind0], x[ind1], alternative="g")$p.value
  })

  pval2 = apply(sci_result[[1]], 1, function(x) {
    if(mean(x[ind0]) < mean(x[ind1])+ 0.05) return(1)
    wilcox.test(x[ind0], x[ind1], alternative="g")$p.value
  })

  plot(log10(res$padj)*(res$log2FoldChange < -1), log10(pval1))

  caTools::colAUC(pval1, DE_t)
  caTools::colAUC(pval2, DE_t)


  pred1 = prediction(pval1, 1-DE_t)
  pred2 = prediction(pval2, 1-DE_t)

  roc.perf1 = performance(pred1, measure = "tpr", x.measure = "fpr")
  roc.perf2 = performance(pred2, measure = "tpr", x.measure = "fpr")

  plot(roc.perf1); plot(roc.perf2, col=2,add=T);
  legend("bottomright",c("MAGIC", "IMPUTE"), col=1:3, border=F,lty=1)


}

gplots::venn(list("faimpute"= which(pval1<0.001), "scimpute" = which(pval2<0.00001), "true" = which(fdr_t<0.05)))

#ggplot(melt(result$beta[1:200,]), aes(Var2,Var1)) +geom_tile(aes(fill = value)) +
#scale_fill_gradient2()+theme_bw()+xlab("")+ylab("")



#### Analyzing results ####
# heatmap for differential expressed genes
source("Heatmap_factormodel.R")
# get differential expressed genes by imputed single cell
# combine H1 and H9 cell types
celltype0 = celltype
celltype = plyr::mapvalues(celltype, from="H9", to="H1")
# wilcox-test using imputed expression
impute_DE <- matrix(0, nrow(Y2), 6)
impute_DE2 <- matrix(0, nrow(Y2), 6)
impute_DE0 <- matrix(0, nrow(Y2), 6)
i=1
for(m in levels(factor(celltype)))
{
  ind0 = which(celltype==m)
  pvs = apply(Y2, 1, function(x) { #sci_result[[1]] ,result2$impt
    if(mean(x[ind0]) < mean(x[-ind0]) + 0.2) return(1)
    wilcox.test(x[ind0], x[-ind0], alternative="g")$p.value #t
  })
  fdr = p.adjust(pvs)
  #impute_DE[fdr<0.001 ,i] = 1
  impute_DE0[order(fdr)[1:100],i] = 1
  i = i+1

}


# heatmap of orignial and imputed dataset
# 3 levels of imputation
imputeY = result$Y*result$posSd + result$posM # Y, impute for 0
imputeY = t(result$Yimp0)*result$posSd + result$posM # E(Y), impute for all data
imputeY = result2$impt # E(Y) impute for 0

imputeY = t(data_MAGIC$result) #

rownames(imputeY) = 1:nrow(imputeY)
colnames(imputeY) = 1:ncol(imputeY)


sel_gene0 = which(rowSums(true_DE) > 0)
sel_gene1 = which(rowSums(impute_DE) > 0)
sel_gene2 = which(rowSums(impute_DE2) > 0)
sel_gene3 = which(rowSums(impute_DE0) > 0)

pdf("chu_data/DE_genes_venn.pdf")
gplots::venn(list("bulk" = sel_gene0, "faimpute" = sel_gene1,  "org" = sel_gene3))#"scimpute" = sel_gene2,
dev.off()

figures2 <- heatmap(imputeY[sel_gene0,], factor(celltype0), factor(celltype), grp_ord = NULL, scale=T,cell.ord ="celltype") # use most imputed data to order gene!, ,
figures <- heatmap(result2$impt[sel_gene0, ], factor(celltype0), factor(celltype), grp_ord = NULL, scale=T,cell.ord ="celltype",gene_ord = figures2[[2]], disp.max = 3)

figures <- heatmap(Y2[sel_gene0, ], factor(celltype0), factor(celltype), grp_ord = NULL, scale=T,cell.ord ="celltype",gene_ord = figures2[[2]], disp.max = 3)

figures2 <- heatmap(sci_result[[1]][sel_gene2,], factor(celltype0), factor(celltype), grp_ord = NULL, scale=T,cell.ord ="celltype") # use most imputed data to order gene!, ,



plot_grid(figures2[[1]][[2]],figures2[[1]][[1]], nrow = 3, align = "v",rel_heights=c(1,7)) #figures[[3]],
pdf("chu_data/scimpute_DE_heatmap.pdf")
plot(figures2[[1]][[1]])
dev.off()
plot(figures[[1]][[1]])



#
g = which(rownames(Y2) == "SOX17")
boxplot(result2$impt[g,]~celltype)
boxplot(sci_result[[1]][g,]~celltype)
boxplot(Y2[g,]~celltype)

# MSE of imputed expression vs. truth
sqrt(mean((result2$impt[Y2==0] - Y[Y2==0])^2))
sqrt(mean((sci_result[[1]][Y2==0] - Y[Y2==0])^2)) #cortex.saver.genes
mean((t(data_MAGIC$result)[Y2==0] - simu_data$Y[Y2==0])^2)

for(m in 1:M0)
{
  mean( rowMeans(result2$impt[, celltype==m]) - mu[,m])
}

plot(mu[,1],rowMeans(result2$impt[, celltype==1]))
plot(mu[,1],rowMeans(sci_result[[1]][, celltype==1]))

# differential expressed genes ROC
library(ROCR)
res0 = matrix(0, G, M0)
res = matrix(0, G, M0)
res2 = matrix(0, G, M0)
res3 = matrix(0, G, M0)
label = matrix(0, G, M0)

for(m in 1:M0)
{
  res0[,m] <- apply(Y2, 1, function(x) { t.test(x~celltype==m, alternative="l")$p.value })
  #res[,m] <- apply(data_MAGIC$result, 2, function(x) { t.test(x~celltype==m, alternative="l")$p.value })
  res2[,m] <- apply(sci_result[[1]], 1, function(x) { t.test(x~celltype==m,alternative="l")$p.value })
  res3[,m] <- apply(result2$impt, 1, function(x) { t.test(x~celltype==m,alternative="l")$p.value })
  # res3[,m] <- apply(1:G, 1, function(i) {
  #   x1 = result2$impt[i, celltype==m]
  #   x2 = result2$impt[i, celltype!=m]
  #   i01 = which(Y2[i, ] ==0 & celltype==m)
  #   i02 = which(Y2[i, ] ==0 & celltype!=m)
  #   i11 = which(Y2[i, ] ==0 & celltype==m)
  #   i12 = which(Y2[i, ] ==0 & celltype!=m)
  #   md = mean(x1) - mean(x2)
  #   sd = sqrt(var(x1[i01])/length(i01) +
  #   })

  ss = ((m-1) * S0 + 1) : (m * S0)
  label[ss ,m] = 1
} #

pred0 = prediction(res0, 1-label)
pred1 = prediction(res, 1-label)
pred2 = prediction(res2, 1-label)
pred3 = prediction(res3, 1-label)
roc.perf0 = performance(pred0, measure = "tpr", x.measure = "fpr")
roc.perf1 = performance(pred1, measure = "tpr", x.measure = "fpr")
roc.perf2 = performance(pred2, measure = "tpr", x.measure = "fpr")
roc.perf3 = performance(pred3, measure = "tpr", x.measure = "fpr")

plot(roc.perf0); plot(roc.perf2, col=2,add=T); plot(roc.perf3, col=3,add=T)
legend("bottomright",c("MAGIC", "SAVER", "IMPUTE"), col=1:3, border=F,lty=1)

for(m in 1:M0) print(caTools::colAUC(res0[,m], label[,m]))
for(m in 1:M0) print(caTools::colAUC(res2[,m], label[,m]))
for(m in 1:M0) print(caTools::colAUC(res3[,m], label[,m]))

# check variance of imputed values
impt_sd = sqrt(result$beta^2 %*% t(result2$varF) + result$sigma[,1])*result$posSd




# correlation
true_cors <- cor(t(simu_data$Y))
impute_cors <- cor(t(result2$impt)) #data_MAGIC$result
impute_cors <- cor(t(cortex.saver.genes)) #data_MAGIC$result
plot(c(true_cors), c(impute_cors));abline(c(0,1))


EF = rep(0, n, 10 )
for(m in 1:ncol(result$z))
{
  EF = EF + result$Ef[[m]] *result$z[,m]
}

pdf("chu_data/faimpute_K6_estimate_factors.pdf")
par(mfrow=c(3,3), mar=c(5,4,1,1), mgp=c(2,0.5,0))
for(k in 1:10) #5
{
  boxplot(EF[,k]*sqrt(mean(result$beta[,k]^2))~celltype2, las=2, outline=F, ylab=paste("Factor",k))
}
dev.off()

K=10
sel_gene = matrix(0, 200, K)
for(k in 1:K)
{
  #print(merge(gene_id[order(abs(result$beta[,k]), decreasing = T)[1:20]], hgnc, by=1))
  sel_gene[,k] = rownames(Y2)[order(abs(result$beta[,k]), decreasing = T)[1:200]]
}

write.table(sel_gene, "chu_data/factor_gene_sym_K6.txt", quote=F,row.names = F, col.names = T,sep="\t")


hist(EF[clus0==7,1])
hist(EF[clus0==3,1])
hist(EF[,2])
