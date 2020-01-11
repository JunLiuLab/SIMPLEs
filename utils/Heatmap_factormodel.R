library(reshape2)
library(dplyr)
library(ggplot2)
library(matrixStats)
library(cowplot)
# order genes for heatmap ng, number of genes per cluster to plot dd, data to
# plot: G*N matrix, pS_thres, posterior prob for selected genes, scale, whether
# to center and standarized each gene if cl0==T, cluster 0 will be put to the
# rightmost of the heatmap dist_method: euclidean or manhattan for cluster group
# means select_gene: if genes selected by p_S is >ng * #clust, select cluster
# specific genes to draw, either by auc or t-test, t-test is slow...  cell.ord:
# order the cell by true cell type
heatmap <- function(dd, celltype, clust, grp_ord = NULL, gene_ord = NULL, scale = F, 
    dist_method = "euclidean", G.legend = "top", cell.ord = "celltype", disp.max = 2.5) {
    
    if (scale) {
        dd = t(scale(t(dd)))
    }
    
    # get mean for each cluster clust = celltype #res$I_hat to remove empty cluster
    y = xtabs(~clust)
    # if(cl0){ y = y[names(y)!=0] # exclude cluster 0 }
    dm = sapply(names(y), function(x) {
        # singleton
        if (sum(clust == x) == 1) {
            dd[, clust == x]
        } else {
            rowMeans(dd[, clust == x])
        }
    })
    
    
    if (is.null(grp_ord)) {
        # reorder each cluster by its mean
        grp_ord <- hclust(dist(t(dm), method = dist_method))
        dend = reorder(as.dendrogram(grp_ord), wts = -y)  #,-6
        grp_ord <- as.hclust(dend)
        
        dm = dm[, grp_ord$order]
        grp_ord = names(y)[grp_ord$order]
        
    } else {
        dm = dm[, grp_ord]
        grp_ord = names(y)[grp_ord]
    }
    
    # cluster and reorder genes
    if (is.null(gene_ord)) {
        gene_ord <- hclust(dist(dd))
        sv = apply(dm, 1, which.max)
        dend = reorder(as.dendrogram(gene_ord), wts = sv)
        gene_ord <- as.hclust(dend)
        gene_ord = gene_ord$order
    }
    
    ## pdf('heatmap_selectmarkers_deng_K9.pdf', height = 9)
    figures <- DoHeatmap(dd[gene_ord, ], clust, celltype, slim.col.label = TRUE, 
        remove.key = F, cex.row = 4, group.order = as.character(grp_ord), rotate.key = T, 
        G.legend = G.legend, cell.ord = cell.ord, disp.max = disp.max)
    return(list(figures, gene_ord))
    ## dev.off()
}



DoHeatmap <- function(data, grp.ident, celltype, disp.min = -2.5, disp.max = 2.5, 
    group.order = NULL, draw.line = TRUE, col.low = "#FF00FF", col.mid = "#000000", 
    col.high = "#FFFF00", slim.col.label = FALSE, remove.key = FALSE, rotate.key = T, 
    title = NULL, cex.col = 10, cex.row = 10, group.label.loc = "bottom", group.label.rot = FALSE, 
    group.cex = 12, group.spacing = 0.15, do.plot = TRUE, G.legend = "top", cell.ord = "celltype") {
    
    PY = colorRampPalette(c("magenta", "yellow"))(50)
    grp.ident = factor(grp.ident, levels = group.order)
    
    # within cluster reorder the cell according to cell type
    celltype <- data.frame(celltype = celltype, ident = grp.ident, group = "0")
    celltype$cell <- colnames(data)
    
    if (cell.ord == "celltype") {
        ord <- order(celltype$ident, celltype$celltype)
    } else if (cell.ord == "nG") {
        ord <- order(celltype$ident, colMeans(data != 0))
    } else {
        ord = 1:nrow(celltype)
    }
    
    celltype <- celltype[ord, ]
    celltype$cell = factor(celltype$cell, levels = celltype$cell)
    
    data = data[, ord]
    
    disp.max = min(disp.max, max(data))
    disp.min = max(disp.min, min(data))
    
    data[data > disp.max] <- disp.max
    data[data < disp.min] <- disp.min
    
    data.use <- as.data.frame(t(data))
    data.use$cell <- rownames(x = data.use)
    colnames(x = data.use) <- make.unique(names = colnames(x = data.use))
    data.use <- data.use %>% melt(id.vars = "cell")
    names(x = data.use)[names(x = data.use) == "variable"] <- "gene"
    names(x = data.use)[names(x = data.use) == "value"] <- "expression"
    data.use$ident <- celltype$ident
    
    # data.use$ident <- factor(data.use$ident, levels = group.order)
    
    breaks <- seq(from = min(data.use$expression), to = max(data.use$expression), 
        length = length(PY) + 1)
    data.use$gene <- with(data = data.use, expr = factor(x = gene, levels = rev(x = unique(x = data.use$gene))))
    data.use$cell <- with(data = data.use, expr = factor(x = cell, levels = colnames(data)))
    if (rotate.key) {
        key.direction <- "horizontal"
        key.title.pos <- "top"
    } else {
        key.direction <- "vertical"
        key.title.pos <- "top"  #'left'
    }
    ### data.use2 <- merge(data.use, celltype[,1:2],by='cell') Expression
    heatmap <- ggplot(data = data.use, mapping = aes(x = cell, y = gene, fill = expression)) + 
        geom_tile() + scale_fill_gradient2(low = col.low, mid = col.mid, high = col.high, 
        name = "", guide = guide_colorbar(direction = key.direction, title.position = key.title.pos)) + 
        scale_y_discrete(position = "right", labels = rev(rownames(data))) + theme(axis.line = element_blank(), 
        axis.title.y = element_blank(), axis.ticks.y = element_blank(), strip.text.x = element_text(size = group.cex), 
        axis.text.y = element_text(size = cex.row), axis.text.x = element_text(size = cex.col), 
        axis.title.x = element_blank())
    if (slim.col.label) {
        heatmap <- heatmap + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), axis.line = element_blank(), axis.title.y = element_blank(), 
            axis.ticks.y = element_blank())
    } else {
        heatmap <- heatmap + theme(axis.text.x = element_text(angle = 90))
    }
    
    if (group.label.loc == "top") {
        switch <- NULL
    } else {
        switch <- "x"
    }
    heatmap <- heatmap + facet_grid(facets = ~ident, drop = TRUE, space = "free", 
        scales = "free", switch = switch) + scale_x_discrete(expand = c(0, 0), drop = TRUE)
    if (draw.line) {
        panel.spacing <- unit(x = group.spacing, units = "lines")
    } else {
        panel.spacing <- unit(x = 0, units = "lines")
    }
    heatmap <- heatmap + theme(strip.background = element_blank(), panel.spacing = panel.spacing)
    if (group.label.rot) {
        heatmap <- heatmap + theme(strip.text.x = element_text(angle = 90))
    }
    
    if (remove.key) {
        heatmap <- heatmap + theme(legend.position = "none")
    } else {
        heatmap <- heatmap + theme(legend.position = "bottom", legend.justification = "center", 
            legend.box.margin = margin(-20, 0, 0, 0))
    }
    if (!is.null(x = title)) {
        heatmap <- heatmap + labs(title = title)
    }
    
    # add number of expressed genes
    nG = data.frame(nG = colMeans(data != 0), cell = celltype$cell, ident = celltype$ident, 
        group = "0")
    mid = median(nG$nG)
    heatmap1 <- ggplot(data = nG, mapping = aes(x = cell, y = group, fill = nG)) + 
        geom_tile() + scale_y_discrete(expand = c(0, 0), position = "right", labels = "") + 
        scale_fill_gradient2(midpoint = mid) + theme(axis.line = element_blank(), 
        axis.title.y = element_blank(), axis.ticks.y = element_blank(), strip.text.x = element_text(size = group.cex), 
        axis.text.y = element_text(size = cex.row), axis.text.x = element_text(size = cex.col), 
        axis.title.x = element_blank()) + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.line = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank()) + facet_grid(facets = ~ident, drop = TRUE, 
        space = "free", scales = "free", switch = switch) + scale_x_discrete(expand = c(0, 
        0), drop = TRUE) + theme(strip.background = element_blank(), panel.spacing = panel.spacing, 
        strip.text.x = element_blank()) + theme(legend.position = G.legend, legend.text = element_text(size = cex.col - 
        1), legend.title = element_text(size = cex.col - 1))
    # add cell type
    heatmap2 <- ggplot(data = celltype, mapping = aes(x = cell, y = group, fill = celltype)) + 
        geom_tile() + scale_y_discrete(expand = c(0, 0), position = "right", labels = "") + 
        theme(axis.line = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), 
            strip.text.x = element_text(size = group.cex), axis.text.y = element_text(size = cex.row), 
            axis.text.x = element_text(size = cex.col), axis.title.x = element_blank()) + 
        theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
            axis.line = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) + 
        facet_grid(facets = ~ident, drop = TRUE, space = "free", scales = "free", 
            switch = switch) + scale_x_discrete(expand = c(0, 0), drop = TRUE) + 
        theme(strip.background = element_blank(), panel.spacing = panel.spacing, 
            strip.text.x = element_blank()) + theme(legend.position = "top", legend.text = element_text(size = cex.col - 
        1), legend.title = element_text(size = cex.col - 1))
    
    if (do.plot) {
        heatmap <- heatmap + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
        heatmap2 <- heatmap2 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
        heatmap1 <- heatmap1 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
        # plot_grid(heatmap2,heatmap, nrow = 2, align = 'v',rel_heights=c(1,7))
    }
    
    return(list(heatmap, heatmap2, heatmap1))
}


# ## test ## load('camp2_seurat.rdat') load('version2/camp2_result_1-12.rdat') #
# plot percentage of expressed genes, order cell within cluster by true celltype
# figures <- heatmap(res$K9, seurat@meta.data$cell_type1)
# plot_grid(figures[[2]],figures[[3]],figures[[1]], nrow = 3, align =
# 'v',rel_heights=c(1,1,7)) # plot percentage of expressed genes, order cell
# within cluster by nG figures <- heatmap(res$K9,
# seurat@meta.data$cell_type1,cell.ord = 'nG')
# plot_grid(figures[[2]],figures[[3]],figures[[1]], nrow = 3, align = 'v',
# rel_heights=c(1,1,7)) # plot percentage of expressed genes (nG), order cell
# within cluster by nG, no legend for nG figures <- heatmap(res$K9,
# seurat@meta.data$cell_type1,cell.ord = 'nG', G.legend = 'none')
# plot_grid(figures[[2]],figures[[3]],figures[[1]], nrow = 3, align =
# 'v',rel_heights=c(1,0.2,7)) heatmap(res, dd=dat@assays$data$logcounts,
# select_gene = 't-test')
