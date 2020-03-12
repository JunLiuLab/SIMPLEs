import::from(here, here)
import::from(ggplot2, ggplot, aes, geom_point, theme_bw, labs, ggtitle)
import::from(Matrix, colSums, rowSums, t)
import::from(matrixStats, rowSds)
import::from(ontologyIndex, get_ontology, get_descendants)
import::from(mclust, adjustedRandIndex)
import::from(Rtsne, Rtsne)

library(ggplot2)
library(RColorBrewer)
library(randomcoloR)
library(ggpubr)
library(gridGraphics)
library(gridExtra)
library(Seurat)
library(cowplot)

tsneplot <- function(tsne_res, types, title="TSNE", palette=NULL,file=NULL, ...) {
  dat <- data.frame(cbind(tsne_res), types)
  colnames(dat) <- c("TSNE1", "TSNE2", "Type")
  if (is.null(palette)){
    palette <- unname(distinctColorPalette(length(unique(types))))
  }

  p <- ggplot(dat, aes(x=TSNE1, y=TSNE2, color=Type)) +
    geom_point() +
    theme_bw()+
    labs(title = title) +
    scale_color_manual(values=palette)
  if (is.null(file)) {
    quartz(title=title, ...)
    plot(p)
  } else {
    pdf(file=file, ...)
    plot(p)
    dev.off()
  }
  return(palette)
}


get_tsne <- function(gene2cell, filename, seed=0, initial_dims = 50) {
  set.seed(0)
  tsne_res <- Rtsne(t(gene2cell), pca_scale=T, pca_center=T, initial_dims=initial_dims,
                    pca=T, check_duplicates=FALSE)
  saveRDS(object=tsne_res$Y, file = filename)
}

load_simple_res <- function(K0=10, M0=2, pm=0.3, cutoff=0.1) {
  prefix <- paste0("SIMPLE_", "K0-",K0, "_M0-", M0,"_pm-", pm)
  simple_res_local_name <- paste0(prefix, "_cufoff-", cutoff, "_maxlambda_result.rds")
  filename <- here(simple_res_local_name)
  if (!file.exists(filename)) {
    stop(paste0(filename, " Not Exists."))
  }
  simple_imput_res <- readRDS(file=filename)
  return(simple_imput_res)
}

simple_get_tsne <- function(K0=10, M0=2, pm=0.3, cutoff=0.1) {
  prefix <- paste0("SIMPLE_", "K0-",K0, "_M0-", M0,"_pm-", pm, "_maxlambda")
  simple_imput_res <- load_simple_res(K0, M0, pm, cutoff)
  set.seed(0)
  simple_tsne_res <- Rtsne(t(simple_imput_res$impt),
                           pca_scale=T, pca_center=T, initial_dims=50,
                           pca=T, check_duplicates=FALSE)
  saveRDS(object = simple_tsne_res$Y,
          file=here(paste0("leukocytes_", prefix,"_tsne.rds")))
}

simple_draw_tsne <- function(meta_data, palettes, K0=10, M0=2, pm=0.3, palette=NULL){
  ## prefix <- paste0("SIMPLE_", "K0-",K0, "_M0-", M0,"_pm-", pm)
  prefix <- paste0("SIMPLE_", "K0-",K0, "_M0-", M0,"_pm-", pm, "_maxlambda")
  local_file_name <- paste0("leukocytes_",prefix, "_tsne.rds")
  filename <- here(local_file_name)
  if (!file.exists(filename)) {
    warning(paste0("No such file: ", filename))
  } else {
    for (typenm in c("tissuecell", "tissue", "cell_ontology_class")){
      width <- 12
      height <- 7
      if (typenm == "tissuecell") {
        width <- 20
        height <- 10
      }
      tsne_res <- readRDS(filename)
      draw_save_file=here(paste("leukocytes", prefix, typenm,"tsne.pdf", sep = "_"))
      types <- get(typenm, meta_data)
      palette <- get(typenm, palettes)
      tsneplot(tsne_res = tsne_res, types = types, title="SIMPLE", palette=palette,
               file = draw_save_file, width=width,height=height)
    }
  }
}

## * raw_data
## Generate Figure 8(b), Figure S11, and Figure S12.
handle_raw_data <- function(){
  normalized_data <- readRDS(here("leukocytes_normalized_data.rds"))
  meta_data <- readRDS(here("leukocytes_meta_data.rds"))
  meta_data$tissuecell <- paste(meta_data$tissue, meta_data$cell_ontology_class)

  ## ** calculate tsne and save
  get_tsne(normalized_data, here("leukocytes_raw_tsne.rds"))

  ## ** load tsne and plot
  raw_tsne_res <- readRDS(file=here("leukocytes_raw_tsne.rds"))

  raw_tc_palette <- tsneplot(raw_tsne_res, types=meta_data$tissuecell, title="raw",
                      file=here("leukocytes_raw_tsne_tc.pdf"),
                      width=20, height=10)
  raw_t_palette <- tsneplot(raw_tsne_res, types=meta_data$tissue, title="raw",
          file=here("leukocytes_raw_tsne_tissue.pdf"),
          width=10,height=7)
  raw_c_palette <- tsneplot(raw_tsne_res, types=meta_data$cell_ontology_class,
                            title="raw",
          file=here("leukocytes_raw_tsne_cell.pdf"),
          width=12, height=7)

  raw_palettes <- list(tissuecell=raw_tc_palette,
                      cell_ontology_class=raw_c_palette,
                      tissue=raw_t_palette)
  # [Important] save the color palettes for later usages.
  saveRDS(object = raw_palettes, file=here("raw_palettes.rds"))
}

## * SIMPLE result analyze
############################## generate the tsne file.
handle_simple_tsne <- function() {
  for (K0 in c(5, 10,15, 20)) {
    for (M0 in c(1, 2, 8, 12, 22)) {
      pm  <- 0.4
      prefix <- paste0("SIMPLE_", "K0-",K0, "_M0-", M0,"_pm-", pm, "_maxlambda")
      tsne_res_file <- here(paste0("leukocytes", prefix,"_tsne.rds"))
      if( file.exists(tsne_res_file)  ) {
        next
      }
      tryCatch({
        message("Current file: ", prefix)
        simple_get_tsne(K0, M0, pm, 0.1)
      },
      error=function(e){
        message(e)
      },
      finally = {
        next
      })
    }
  }
}

## draw the scatter plot based tsne result.
############################### generate Figure 8a.
handle_simple_draw_tsne <- function() {
  meta_data <- readRDS(here("leukocytes_meta_data.rds"))
  meta_data$tissuecell <- paste(meta_data$tissue, meta_data$cell_ontology_class)
  palettes <- readRDS(here("raw_palettes.rds"))
  # do tsne
  for (K0 in c(5, 10, 15, 20)) {
    for (M0 in c(1, 2, 8,12,22)) {
      pm <- 0.4
      simple_draw_tsne(meta_data=meta_data, palettes=palettes,K0=K0, M0=M0, pm=pm)
    }
  }
}


## * B metrix analysis
### draw the gene mdoules according to different latent factors.
meta_data <- readRDS(here("leukocytes_meta_data.rds"))
meta_data$tissuecell <- paste(meta_data$tissue, meta_data$cell_ontology_class)
palettes <- readRDS(here("raw_palettes.rds"))
######################## generate Figure 9(a), and Figure S13.
draw_simple_boxplot <- function(K0=10, M0=12, pm=0.4, cutoff=0.1,
                                palettes=palettes, onwhat="cell",
                                oncondition=NULL, condname ="allcell",
                                width = 14, height=40, outlier.shape=NA
                                ) {
  simple_res <- load_simple_res(K0, M0, pm, cutoff)
  ## sort(colMeans(abs(simple_res$beta)))
  B <- simple_res$beta
  Ef <- data.frame(simple_res$Ef)
  colnames(Ef) <- paste("Factor", c(1:K0), sep = ".", collapse = NULL)
  ## colnames(B) <- paste("Factor", c(1:K0), sep = ".", collapse = NULL)
  for (k in seq(1,K0)) {
    Ef[, k] <- Ef[, k] * sqrt(mean(B[, k]^2))
  }
  Ef$cell <- factor(meta_data$cell_ontology_class)
  Ef$tissue <- factor(meta_data$tissue)

  if(!is.null(oncondition)){
    if (onwhat == "cell") {
      Ef <- Ef[Ef$tissue == oncondition, ]
      print(dim(Ef))
    }
    if (onwhat == "tissue") {
      Ef <- Ef[ is.element(Ef$cell, oncondition), ]
      print(dim(Ef))
    }
  }
  ## sort by the median of tissue.
  pltList <- lapply(1:K0, function(k) {
    p <- ggplot(Ef, aes_string(
                      x= paste0("reorder(",onwhat,",-",paste("Factor", k, sep="."),
                               ",FUN = median)"),
                              y=paste("Factor", k, sep = "."),
                              fill=onwhat
                              )) +
      ## ylim(-1,1) +
      geom_boxplot(outlier.shape = outlier.shape) +
      theme_classic() +
      scale_fill_manual(values=palettes) +
      theme(
        axis.text.x=element_text(angle=45, color="black", size = 16, family="Helvetica"),
        axis.text.y=element_text(family = "Helvetica", color="black", size=10),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=16, family="Helvetica"),
        )  +
      ylab(paste("FACTOR", k, sep = " "))
  } )

  prefix <- paste0("SIMPLE_", "K0-",K0, "_M0-", M0,"_pm-", pm)
  if (!is.null(oncondition)) {
    prefix <- paste0(prefix, "_", condname)
  }
  boxplot_save_fnm <- here(paste0(prefix, "_",onwhat, "_boxplot.pdf"))
  ggexport(plotlist = pltList, filename = boxplot_save_fnm,
          nrow = 10, ncol = 1, width = width, height=height)
}


## ** modify seurat package functions and do dotplots.
PercentAbove <- function(x, threshold) {
  return(length(x = x[x > threshold]) / length(x = x))
}

`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

MyDotPlot <- function(object, assay = NULL, features, cols = c("lightgrey","blue"),
                      col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
                      group.by = NULL, split.by = NULL, scale.by = "radius",
                      scale.min = NA, scale.max = NA) {
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    scale.func <- switch(EXPR = scale.by, size = scale_size, 
        radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
    data.features <- FetchData(object = object, vars = features)
    data.features$id <- if (is.null(x = group.by)) {
        Idents(object = object)
    }
    else {
        object[[group.by, drop = TRUE]]
    }
    if (!is.factor(x = data.features$id)) {
        data.features$id <- factor(x = data.features$id)
    }
    id.levels <- levels(x = data.features$id)
    data.features$id <- as.vector(x = data.features$id)
    if (!is.null(x = split.by)) {
        splits <- object[[split.by, drop = TRUE]]
        if (length(x = unique(x = splits)) > length(x = cols)) {
            stop("Not enought colors for the number of groups")
        }
        cols <- cols[1:length(x = unique(x = splits))]
        names(x = cols) <- unique(x = splits)
        data.features$id <- paste(data.features$id, splits, sep = "_")
        unique.splits <- unique(x = splits)
        id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
            "_", rep(x = unique(x = splits), times = length(x = id.levels)))
    }
    data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
        data.use <- data.features[data.features$id == ident, 
            1:(ncol(x = data.features) - 1), drop = FALSE]
        avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
            ## return(mean(x = expm1(x = x)))
            return(mean(x = x))
        })
        pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
            threshold = 0)
        return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    })
    names(x = data.plot) <- unique(x = data.features$id)
    data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
        data.use <- as.data.frame(x = data.plot[[x]])
        data.use$features.plot <- rownames(x = data.use)
        data.use$id <- x
        return(data.use)
    })
    data.plot <- do.call(what = "rbind", args = data.plot)
    if (!is.null(x = id.levels)) {
        data.plot$id <- factor(x = data.plot$id, levels = id.levels)
    }
    avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
        FUN = function(x) {
            data.use <- data.plot[data.plot$features.plot == x, "avg.exp"]
            ## * Change to no scale
            ## data.use <- scale(x = data.use)
            data.use <- MinMax(data = data.use, min = col.min, max = col.max)
            return(data.use)
        })
    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    if (!is.null(x = split.by)) {
        avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
            breaks = 20))
    }
    data.plot$avg.exp.scaled <- avg.exp.scaled
    data.plot$features.plot <- factor(x = data.plot$features.plot, 
        levels = rev(x = features))
    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (!is.null(x = split.by)) {
        splits.use <- vapply(X = strsplit(x = as.character(x = data.plot$id), 
            split = "_"), FUN = "[[", FUN.VALUE = character(length = 1L), 
            2)
        data.plot$colors <- mapply(FUN = function(color, value) {
            return(colorRampPalette(colors = c("grey", color))(20)[value])
        }, color = cols[splits.use], value = avg.exp.scaled)
    }
    color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled", 
        no = "colors")
    if (!is.na(x = scale.min)) {
        data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }
    plot <- ggplot(data = data.plot,
                   mapping = aes_string(x = "features.plot",y = "id")) +
      geom_point(mapping = aes_string(size = "pct.exp",color = color.by)) +
      scale.func(range = c(0, dot.scale),limits = c(scale.min, scale.max)) +
      theme(axis.title.x = element_blank(),axis.title.y = element_blank()) +
      guides(size = guide_legend(title = "Percent Expressed")) +
      labs(x = "Features", y = ifelse(test = is.null(x = split.by), 
          yes = "Identity", no = "Split Identity")) + theme_cowplot()
    if (!is.null(x = split.by)) {
        plot <- plot + scale_color_identity()
    }
    else if (length(x = cols) == 1) {
        plot <- plot + scale_color_distiller(palette = cols)
    }
    else {
      plot <- plot + scale_color_gradient(low = cols[1], high = cols[2],
                                          breaks=c(0, 5, col.max),
                                          labels=c(0, 5, col.max),
                                          limits=c(0, col.max)
                                          )
    }

    if (is.null(x = split.by)) {
      plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
      ## plot <- plot + scale_fill_gradient(breaks=c(0,col.max, 0.2))
    }
    return(plot)
}

MyPubDotPlot <- function(object, features,
                         group.by,cols,
                         col.max , col.min ,
                         scale.min, scale.max, dot.scale){
  p <- MyDotPlot(object = object, features = features,
                 group.by = group.by,cols = cols,
                 col.max = col.max, col.min = col.min,
                 scale.min = scale.min, scale.max = scale.max)
  p <- p +
    theme(axis.text.x=element_text(angle=90, color="black",
                                   size = 15, family="Helvetica"),
          axis.text.y=element_text(family = "Helvetica", color="black", size=15,
                                   angle = 30
                                   ),
          axis.ticks.x =  element_blank(),
          axis.title = element_blank(),
          axis.line.x =  element_blank()
          )
  return(p)
}


################## generate Figure 9(b)
normalized_data <- readRDS(here("TabulaMuris","data","leukocytes_normalized_data.rds"))
meta_data <- readRDS(here("TabulaMuris","data","leukocytes_meta_data.rds"))
meta_data$tissuecell <- paste(meta_data$tissue, meta_data$cell_ontology_class)

TcellFamily <- c("T cell", "immature T cell","regulatory T cell",
                 "immature NK T cell")
BcellFamily <- c("immature B cell", "B cell", "naive B cell",
                 "precursor B cell","late pro-B cell")

## * stats cells in tissue
tcell_index <- is.element(meta_data$cell_ontology_class, TcellFamily)
tcells <- meta_data$cell_ontology_class[tcell_index]
tcell_tissues <- meta_data$tissue[tcell_index]
table(tcells, tcell_tissues)

## * choose the specific parameters
M0 <- 1
K0 <- 10
pm <- 0.4
cutoff <- 0.1

## * load simple result
simple_res <- load_simple_res(K0, M0, pm, cutoff)
B <- simple_res$beta
colnames(B) <- paste("Factor", c(1:K0), sep = ".", collapse = NULL)
rownames(B) <- rownames(normalized_data)

## * construct seurat object
raw_seurat <- CreateSeuratObject(
  counts = normalized_data[,tcell_index],
  project = "RAW",assay = "imputation")
raw_seurat <- AddMetaData(object = raw_seurat, metadata = tcell_tissues,
                             col.name = "tissue")
raw_seurat <- AddMetaData(object = raw_seurat,
                             metadata = tcells,
                          col.name = "cell")

simple_seurat <- CreateSeuratObject(
  counts = simple_res$impt[,tcell_index],
  project = "SIMPLE",assay = "imputation")
simple_seurat <- AddMetaData(object = simple_seurat, metadata = tcell_tissues,
                             col.name = "tissue")
simple_seurat <- AddMetaData(object = simple_seurat,
                             metadata = tcells,
                             col.name = "cell")

get_top_genes <- function(factor = 9, ph = 10) {
  factornm <- paste("Factor", factor, sep=".")
  gene_module <- B[, colnames(B) == factornm]
  Bbi <- order(gene_module, decreasing = T)
  Bb <- gene_module[Bbi]
  top_high_genes <- names(Bb)[1:ph]
  return(top_high_genes)
}

################## generate Figure 9(b)
f3_top_genes <- get_top_genes(factor=3)
f9_top_genes <- get_top_genes(factor=9)

scale.min <- 15
scale.max <- 100
col.max <- 10


orig_p_f3 <- MyPubDotPlot(object=raw_seurat, features=f3_top_genes,
                       group.by = "tissue", cols = c("lightgrey", "blue"),
                       col.max = col.max, col.min = 0,
                       scale.min = scale.min,scale.max = scale.max)
simp_p_f3 <- MyPubDotPlot(object = simple_seurat, features = f3_top_genes,
                    group.by = "tissue",cols = c("lightgrey","blue"),
                    col.max = col.max, col.min = 0,
                    scale.min = scale.min,scale.max = scale.max)

orig_p_f9 <- MyPubDotPlot(object=raw_seurat, features=f9_top_genes,
                          group.by = "tissue", cols = c("lightgrey", "blue"),
                          col.max = col.max, col.min = 0,
                          scale.min = scale.min,scale.max = scale.max)

simp_p_f9 <- MyPubDotPlot(object = simple_seurat, features = f9_top_genes,
                          group.by = "tissue",cols = c("lightgrey","blue"),
                          col.max = col.max, col.min = 0,
                          scale.min = scale.min,scale.max = scale.max)

os_p <- ggarrange(orig_p_f3, simp_p_f3, orig_p_f9, simp_p_f9,
                  nrow=2, ncol=2, legend = "right", 
                  common.legend = TRUE)
os_p
## os_p <- ggarrange(orig_p, simp_p, nrow=1, ncol=2, legend = "none", widths = c(1, 1))
f_figure_nm <- paste0("f-", f, "_dotplot.pdf")
pdf(file = here(f_figure_nm),
    width = 8,height = 4)
print(os_p)
dev.off()
