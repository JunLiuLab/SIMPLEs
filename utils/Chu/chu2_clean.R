library(scImpute)
library(Rmagic)
library(foreach)
library(doParallel)
library(SIMPLE)
library(matrixStats)
source("Simulation/utils.R")
load("../../chu_data/data_ts.rdat")

#### imputation ####
K0=10
M0=6
pm=0.8
n = ncol(dat_norm)
Y2 = dat_norm
sel = which(rowSds(Y2)>1.5)
Y2 = Y2[sel, ]
bulk_mean = bulk_norm_mean[sel,]
celltype_true = celltype

#DE = DE[sel,]

getCluster(Y2, celltype, Ks=c(5, 10,15, 20))
Y2_scale = t(scale(t(Y2)))
s = svd(Y2_scale)
km0 <- kmeans(t(Y2_scale)%*% s$u[,1:15], M0, iter.max = 80, nstart = 300) #s$v[,1:K]
#km0 <- kmeans(s$v[,1:20], M0, iter.max = 80, nstart = 300)
clus = km0$cluster

registerDoParallel(cores = 7)
result <- scimpclu(Y2, K0, M0, clus = clus, K = 15, iter= 10, est_z = 10, impt_it = 1, max_lambda=T, est_lam = 1, penl = 1, sigma0 = 100, p_min = pm, cutoff=0.01,verbose=T, min_gene = 500, num_mc=3) #as.numeric(factor(celltype_true)) NULL, est_z = 10 clus
resultY = result$Y

result1 = do_impute(Y2, result$Y, result$beta, result$lambda, result$sigma, result$mu, result$pi, result$geneM, result$geneSd, result$initclus, mcmc=50, burnin = 5, verbose=F, pg = result$pg, cutoff = 0.01) # only for AUC and clustering

resultB <- scimpclu_bulk(Y2, K0, bulk_mean, M0, celltype=as.numeric(factor(celltype)), clus = clus, K = 15, iter= 10, est_z = 10, impt_it = 1, max_lambda=T, est_lam = 2, penl =1, sigma0 = 100, p_min = pm, min_gene = 500,cutoff=0.01,verbose=T, num_mc=3) # clus
resultYB = resultB$Y
xtabs(~apply(resultB$z, 1,which.max) + celltype)

result2 = do_impute(Y2, resultB$Y, resultB$beta, resultB$lambda, resultB$sigma, resultB$mu, resultB$pi, resultB$geneM, resultB$geneSd, rep(1, n), mcmc=50, burnin = 5, verbose=F, pg = resultB$pg, cutoff = 0.01)

getCluster(result2$impt, celltype_true, Ks = c(10,15,20), M0 = M0)

# result <- doEM_impute_bulk(Y2, 6, K0, bulk_mean, as.numeric(factor(celltype)), iter=5, est_lam = 3, verbose=T, max_lambda=F, penl =2, clus = NULL, est_z = 1,sigma0 = 100, p_min =0.7,cutoff = 0.1)
#
# km <- kmeans(result$Ef[[1]], 6, iter.max = 80, nstart = 300)
# xtabs(~km$cluster + celltype)

#result2 = do_impute(Y2, result$Y, result$geneM, result$geneSd, 1, K0, result$beta, result$lambda, result$sigma, result$mu, result$pi, mcmc=90, burnin = 10, verbose=F, cutoff = 0.1, pg=result$pg, clus = as.numeric(factor(celltype)))

#save(result, result2, file="factor_model_result/chu_ts_bulk_0228-099.rdat")
save(result1, result,file="../factor_model_result/chu_ts_SIMPLE_0710.rdat")
save(result1, result2, result, resultB, file="factor_model_result/chu_ts_SIMPLE_bulk_0710.rdat")


# not include bulk
result <- doEM_impute(Y2, 6, K0, celltype = celltype, iter=5, est_lam = 3, verbose=T, max_lambda=F, penl =2, clus = NULL, est_z = 1,sigma0 = 100, p_min =0.7, cutoff = 0.1, K = 20)
result2 = do_impute(Y2, result$Y, result$geneM, result$geneSd, 1, K0, result$beta, result$lambda, result$sigma, result$mu, result$pi, pg = result$pg, mcmc=90, burnin = 10, verbose=T, cutoff = 0.1)


sci_result <- scimpute("", infile = "csv", outfile = "csv", type = "count",
                       "../scimpute/", labeled = FALSE, drop_thre = pm, Kcluster = 6,
                       labels = NULL, genelen = NULL, ncores =6, count_lnorm = Y2+log10(1.001)) #0.7

save(sci_result, file="../../chu_data/scimpute_ts_K6_0710.rdat")

MAGIC <- magic(t(Y2), genes="all_genes", t="auto")

getCluster(sci_result[[1]], celltype, Ks=c(10,15, 20))

#tsne

tsneplot(t(Y2), celltype, scale=T, center=T, K=20, file="chu_st_org_tsne.pdf")
tsneplot(t(result2$impt), celltype, scale=T, center=T, K=15, file="chu_ts_myimpute_tsne.pdf") #impute_bulk

tsneplot(t(result1$impt), celltype, scale=T, center=T, K=15, file="chu_ts_myimpute_tsne.pdf")

tsneplot(t(sci_result[[1]]), celltype, scale=T, center=T, K=20, file="chu_ts_scimpute_tsne.pdf")



# correlation with bulk
# correlation between sc and bulk
corB = rep(0,5)
i = 1
for( x in levels(factor(celltype)))
{
  print(x)
  corB[4] = corB[4] + cor(bulk_mean[,i], rowMeans(sci_result[[1]][,celltype == x]))
  #corB[3] = corB[3] + cor(bulk_mean[,i], colMeans(MAGIC$result[celltype == x,]))
  corB[5] = corB[5] + cor(bulk_mean[,i], rowMeans(result2$impt[, celltype == x]))
  corB[2] = corB[2] + cor(bulk_mean[,i], rowMeans(Y2[, celltype == x]))
  i = i+1
}

sci_mean <- sapply(levels(factor(celltype)),function(x) rowMeans(sci_result[[1]][,celltype == x]))
my_mean <- sapply(levels(factor(celltype)),function(x) rowMeans(result2$impt[,celltype == x]))
sci_cor <- sapply(1:nrow(bulk_mean), function(g) cor(bulk_mean[g,], sci_mean[g,]))
my_cor <- sapply(1:nrow(bulk_mean), function(g) cor(bulk_mean[g,],my_mean[g,]))

plot(density(my_cor, na.rm=T))
lines(density(sci_cor, na.rm=T),col=2)
# heatmap for differential expressed genes
source("Heatmap_factormodel.R")
# get differential expressed genes by imputed single cell


# wilcox-test using imputed expression
# impute_DE <- matrix(0, nrow(Y2), 6)
# impute_DE2 <- matrix(0, nrow(Y2), 6)
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
#imputeY = result$Y*result$posSd + result$posM # Y, impute for 0
imputeY = t(result$Yimp0)*result$geneSd + result$geneM # E(Y), impute for all data
#imputeY = result2$impt # E(Y) impute for 0

imputeY = t(data_MAGIC$result) #

rownames(imputeY) = 1:nrow(imputeY)
colnames(imputeY) = 1:ncol(imputeY)


# sel_gene0 = which(rowSums(true_DE) > 0)
# sel_gene1 = which(rowSums(impute_DE) > 0)
# sel_gene2 = which(rowSums(impute_DE2) > 0)
sel_gene3 = which(rowSums(impute_DE0) > 0)

pdf("chu_data/DE_genes_venn.pdf")
gplots::venn(list("bulk" = sel_gene0, "faimpute" = sel_gene1,  "org" = sel_gene3))#"scimpute" = sel_gene2,
dev.off()

figures2 <- heatmap(imputeY[sel_gene3,], factor(celltype), factor(celltype), grp_ord = NULL, scale=T,cell.ord ="celltype") # use most imputed data to order gene!, ,
figures <- heatmap(result2$impt[sel_gene3, ], factor(celltype), factor(celltype), grp_ord = NULL, scale=T,cell.ord ="celltype",gene_ord = figures2[[2]], disp.max = 3)

figures <- heatmap(Y2[sel_gene3, ], factor(celltype), factor(celltype), grp_ord = 1:6, scale=T,cell.ord ="celltype",gene_ord = figures2[[2]], disp.max = 3)

figures2 <- heatmap(sci_result[[1]][sel_gene2,], factor(celltype0), factor(celltype), grp_ord = NULL, scale=T,cell.ord ="celltype") # use most imputed data to order gene!, ,



 #figures[[3]],
#pdf("chu_data/scimpute_DE_heatmap.pdf")
png("factor_model_result/chu_ts_heatmap.png")
plot(figures[[1]][[1]])

plot_grid(figures2[[1]][[2]],figures2[[1]][[1]], nrow = 3, align = "v",rel_heights=c(1,7))
dev.off()

# monocle
library(monocle)
Y2_scale = t(scale(t(Y2)))
monocle_org <- DDRTree(Y2_scale)

metacell <- data.frame("sample"=colnames(Y2),
                       "celltype" = celltype,
                       row.names = colnames(Y2))


pd <- new("AnnotatedDataFrame", data = metacell)
PS2_Mon <- newCellDataSet(as.matrix(Y2), phenoData = pd, featureData = NULL, expressionFamily=VGAM::tobit(),lowerDetectionLimit= -Inf) #sci_result[[1]]

PS2_Mon@expressionFamily@vfamily = "Tobit"
PS2_Mon <- reduceDimension(PS2_Mon, method = 'DDRTree',  verbose = T, norm="none", pseudo_expr = 0)

GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$celltype)[,"12h"]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  }else {
    return (1)
  }
}

PS2_Mon<- orderCells(PS2_Mon)
PS2_Mon <- orderCells(PS2_Mon, root_state = GM_state(PS2_Mon))

plot_cell_trajectory(PS2_Mon, cell_size = 1, show_branch_points =T, color_by = "celltype")
saveRDS(PS2_Mon, file="../chu_data/chu_ts_monocle.rds")
# get genes for each cluster
#cellorder <- pData(PS2_Mon)[ with(pData(PS2_Mon), order(celltype,Pseudotime)),]
temp = pData(PS2_Mon)
temp[temp$celltype=="96h","celltype"]="72h"
ord =  with(temp, order(celltype,Pseudotime))
#cols = c(rgb(0,0,1,0.5), rgb(0,1,0,0.5), rgb(1,0,0,0.5))

# marker genes:
marker=list("12h" = c("NODAL", "EOMES", "ID1"), "24h" = c("T", "MSX2", "CDX1"), "36h" = c("CER1", "GATA4"), "72h" = c("DKK4", "MYCT1"), "96h" = c("PRDM1", "POU2AF1"))

bord = cumsum(xtabs(~celltype))[1:4]
par(mfrow=c(2,2), mar=c(2,2,2,1), cex=1.2, cex.lab=0.9, tcl = -0.2, mgp=c(1,0.2,0))
for(g in unlist(marker))
{
  # check zero count
  ##ord = order(pData(PS2_Mon)$Pseudotime)
  #ind = which(Y2[g,ord ] ==0)
  if(!(g %in% rownames(Y2))) next
  col = (Y2[g,ord ] ==0) +3
  y2 = Y2[g, ord]
  imp = result2$impt[g, ord]
  #imp = sci_result[[1]][g, ord]

  pdf(paste0("factor_model_result/chu_ts_marker_",g,".pdf"))

  scatter.smooth(y2, lpars =list(lwd=2), cex=0.5, pch=16, col=col, xlab="Cells", ylab="Normalized Log UMI", main=g)
  abline(v=bord, lwd=2, lty=2,col="gray")


  scatter.smooth(imp, lpars =list(lwd=2), cex=0.5, pch=16, col=col, xlab="Cells", ylab="Normalized Log UMI", main="Our Method") #g
  abline(v=bord, lwd=2, lty=2,col="gray")
  #scatter.smooth(imp, lpars =list(lwd=2), cex=0.5, pch=16, col=cellorder$celltype, xlab="Cells", ylab="Normalized Log UMI", main=g)

  #scatter.smooth(y2[-ind], lpars =list(lwd=2), cex=0.5, pch=16, col=cellorder$celltype[-ind], xlab="Cells", ylab="Normalized Log UMI", main=g)


  #scatter.smooth(sci_result[[1]][g, ord[ind]], lpars =list(lwd=2), cex=0.5, pch=16, col=cellorder$celltype[ind], xlab="Cells", ylab="Normalized Log UMI", main=g)

  scatter.smooth(sci_result[[1]][g, ord], lpars =list(lwd=2), cex=0.5, pch=16, col=col, xlab="Cells", ylab="Normalized Log UMI", main="scimpute")
  abline(v=bord, lwd=2, lty=2,col="gray")
  scatter.smooth(t(MAGIC$result)[g, ord], lpars =list(lwd=2), cex=0.5, pch=16, col=col, xlab="Cells", ylab="Normalized Log UMI", main="MAGIC")

  abline(v=bord, lwd=2, lty=2,col="gray")

  dev.off()
}


pdf("factor_model_result/chu_ts_gene_corr.pdf")
par(mfrow=c(2,2), mar=c(2,2,2,1), cex=1.2, cex.lab=0.9, tcl = -0.2, mgp=c(1,0.2,0))
#rownames(result2$impt) = rownames(Y2)
plot(Y2["CER1",],Y2["GATA4",], xlab="CER1", ylab="GATA4",main="Data",col=rgb(0,0,1,0.5),pch=16)
legend(-1.5,6.5, paste("corr=",round(cor(Y2["CER1",],Y2["GATA4",]),2)), bty = 'n')
plot(sci_result[[1]]["CER1",], sci_result[[1]]["GATA4",],xlab="CER1", ylab="GATA4",main="scimpute",col=rgb(1,0,0,0.5),pch=16)
legend(-1.5,6.5, paste("corr=",round(cor(sci_result[[1]]["CER1",], sci_result[[1]]["GATA4",]),2)), bty = 'n')
plot(result2$impt["CER1",], result2$impt["GATA4",],xlab="CER1", ylab="GATA4",main="our method",col=rgb(0,1,0,0.5),pch=16)
legend(-3.5,6.5, paste("corr=",round(cor(result2$impt["CER1",], result2$impt["GATA4",]),2)), bty = 'n')
plot(t(MAGIC$result)["CER1",], t(MAGIC$result)["GATA4",],xlab="CER1", ylab="GATA4",main="MAGIC",col=rgb(0,1,1,0.5),pch=16)
legend(-1,4, paste("corr=",round(cor(t(MAGIC$result)["CER1",], t(MAGIC$result)["GATA4",]),2)), bty = 'n')
#plot(Y2["NODAL",], Y2["EOMES",])

dev.off()
