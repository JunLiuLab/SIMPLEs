library(here)
library(ggplot2)
library(matrixStats)
library(tidyverse)

## * configs
others <- c("control", "saver", "rmagic", "scvi", "scrabble", "viper", "scimpute")
pcs <- seq(2, 50, 2)
## * load data
load_de_result <- function(modelnm, k = 2, m = 1, t = "t", thres = "2.0", ...) {
  fnm <- stringr::str_c(modelnm, k, m, t, thres, ".rds", sep = "_")
  readRDS(here::here("10xGenomics", "result", fnm))
}

load_cl_result <- function(modelnm, k = 2, m = 1, ...) {
  t <- "wilcox"
  fnm <- stringr::str_c(modelnm, k, m, t, ".rds", sep = "_")
  r <- readRDS(here::here("10xGenomics", "result", fnm))
  return(r$cl_ari_umi)
}

load_all_cl <- function(...) {
  pcs <- seq(2, 50, 2)
  cl_others <- others %>%
    map(function(x) {
      load_cl_result(x, 2, 1)
    }) %>%
    do.call(rbind, args = .)

  cl_simple <- load_cl_result("simple", 10, 1)
  cls <- rbind(cl_others, cl_simple)
  rownames(cls) <- c(others, "simple")
  colnames(cls) <- paste0("c", pcs)
  print(cls)
  return(cls)
}


## * for differential gene analysis

sum_all_de <- function(t = "t", thres = "2.0") {
  ## * for other methods
  de_others <- others %>%
    map(function(x) {
      load_de_result(x, 2, 1, t, thres)
    }) %>%
    do.call(rbind, args = .)
  de_simple <- load_de_result("simple", 10, 1, t, thres)
  des <- rbind(de_others, de_simple)
  rownames(des) <- c(others, "simple")
  print(des)

  mysum <- data.frame(
    model = rownames(des), mean = rowMeans(des),
    std = rowSds(des)
  )
  print(mysum)
  p <- mysum %>%
    ggplot(data = ., aes(x = model, y = mean, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_errorbar(aes(ymin = mean - std, ymax = mean + std),
      width = .2, position = position_dodge(.9)
    ) +
    theme(
      axis.text.x = element_text(
        color = "black", size = 20,
        family = "Helvetica", face = "bold"
      ),
      axis.text.y = element_text(color = "black", family = "Helvetica", size = 13),
      axis.title.x = element_blank(),
      axis.title.y = element_text(color = "black", family = "Helvetica", size = 20),
      title = element_text(
        family = "Helvetica", color = "black", size = 24,
        face = "bold"
      ),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    ) +
    ylab("AUC")
  prefix <- "AUC of differential expressed gene analysis on \n the dataset sc_10x_5cl with"
  p <- p +
    ggtitle(str_c(prefix, t, "test", sep = " "))
  return(p)
}

a_t2 <- sum_all_de(t = "t", thres = 2.0)
a_w2 <- sum_all_de(t = "wilcox", thres = 2.0)
ggsave(
  path = here::here("10xGenomics", "pdfs"), filename = "de_t2_all.pdf",
  device = "pdf", plot = a_t2, width = 11, height = 11
)
ggsave(
  path = here::here("10xGenomics", "pdfs"), filename = "de_w2_all.pdf",
  device = "pdf", plot = a_w2, width = 11, height = 11
)

## * for clustering analysis
clt_all <- load_all_cl()
clt_all <- tibble::rownames_to_column(as.data.frame(clt_all), "methods")
rownames(clt_all) <- clt_all$methods
p <- ggplot(data = clt_all, aes(
  x = methods, y = c10,
  fill = methods
)) +
  geom_bar(stat = "identity") +
  theme(
    axis.text.x = element_text(
      color = "black", size = 15, family = "Helvetica",
      angle = 90
    ),
    axis.text.y = element_text(color = "black", family = "Helvetica", size = 13),
    axis.title.x = element_blank(),
    axis.title.y = element_text(color = "black", family = "Helvetica", size = 20),
    title = element_text(
      family = "Helvetica", color = "black", size = 24,
      face = "bold"
    ),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  ylab("aRI") +
  ggtitle("clustering analysis on \n the dataset sc_10x_5cl")

ggsave(
  path = here::here("10xGenomics", "pdfs"), filename = "cl_all.pdf",
  device = "pdf", plot = p, width = 11, height = 11
)
