
import::from(here, here)
import::from(Matrix, colSums, rowSums, t)
import::from(matrixStats, rowSds)
import::from(ontologyIndex, get_ontology, get_descendants)

## * load TM counts and metadata.
raw_data <- readRDS(file= here::here("TM_facs_mat.rds" ))
meta_data <- read.csv(
  file=here::here("TM_facs_metadata.csv"),
  row.names = 1,stringsAsFactors = FALSE, header = TRUE, sep = ",")

## ** data preprocessing Remove ERCC and MT genes.
erccs <- grep(pattern = "^ERCC-", x = rownames(x = raw_data), value = TRUE)
percent_ercc <- colSums(raw_data[erccs, ])/colSums(raw_data)
ercc_index <- grep(pattern = "^ERCC-", x = rownames(x = raw_data), value = FALSE)

mts <- grep(pattern = "^MT-", x = rownames(x = raw_data), value = TRUE)
percent_mt <- colSums(raw_data[mts, ])/colSums(raw_data)
mts_index <- grep(pattern = "^MT-", rownames(raw_data), value = FALSE)

raw_data <- raw_data[-c(ercc_index, mts_index), ]

## ** filter low/high  total counts or number of genes.
total_RNA_count <- colSums(raw_data)
count_genes <- colSums(raw_data > 0)

choose_index <- which( (count_genes > 500) & (count_genes < 20000)
                      & (total_RNA_count > 50000) & (total_RNA_count < 2000000) )

facs_raw_data <- raw_data[, choose_index]
facs_meta_data <- meta_data[choose_index, ]

## * select specific cells and data normalization.
get_normalized_data <- function(raw_data_filter_cell, meta_data_filter_cell, ids, 
    selected_genes = NULL,low_percent_count = 0.05, high_percent_count = 0.6, least_std = 1.5, 
    scale_level = 1e+06, margin = 1, print_before_process = FALSE) {
    tcell_index = meta_data_filter_cell$cell_ontology_id %in% ids
    if (print_before_process) {
        print("tissue stats before chosen by gene ontology")
        print(table(factor(meta_data_filter_cell$tissue)))
        print("cell classed stats before chosen by gene ontology")
        print(table(factor(meta_data_filter_cell$cell_ontology_class)))
    }
    tcell_meta <- meta_data_filter_cell[tcell_index, ]

    print("tissue stats after chosen by gene ontoglogy")
    print(table(factor(tcell_meta$tissue)))
    print("cell classed stats after chosen by gene ontology")
    print(table(factor(tcell_meta$cell_ontology_class)))

    num_tcell <- length(which(tcell_index == TRUE))
    pre_tcell_matrix <- raw_data_filter_cell[, tcell_index]
    num_countofcells <- rowSums(pre_tcell_matrix > 0)
    tcell_matrix <- pre_tcell_matrix[which(
        (num_countofcells > low_percent_count * num_tcell)
        & (num_countofcells < high_percent_count * num_tcell) ), ]
    tcell_normalized_data <- log(scale(tcell_matrix, center = FALSE,
                                       scale = colSums(tcell_matrix)) * scale_level + margin)
    par(mfrow=c(2,1))
    boxplot(num_countofcells, main = "boxplot of number of counts of cells genes express")
    boxplot(rowSds(tcell_normalized_data), main="boxplot of standard variance for scaled locally expression")

    final_data_matrix <- as.matrix(tcell_normalized_data[
      which(rowSds(tcell_normalized_data) > least_std), ])

    print(dim(final_data_matrix))

    if (!is.null(selected_genes)) {
      fished_genes  <- grep(paste(selected_genes, collapse = "|"), rownames(final_data_matrix), value = TRUE)
      print(sprintf("fish %d genes from the %d selected genes", length(fished_genes), length(selected_genes)))
      print(fished_genes)
    }

    count_data  <- exp(as.matrix(final_data_matrix)) - 1.0
    count_data[count_data == 0] <- 0.0

    list(normalizded_data = final_data_matrix, meta_data = tcell_meta,
         count_data = count_data)
}

## * All the immune cells
leukocyte_id <- "CL:0000738"
cell_ontologies <- get_ontology(file = "cell_ontology.obo", extract_tags = "everything")
leukocytes <- get_descendants(ontology=cell_ontologies, roots = leukocyte_id,
                              exclude_roots = FALSE)
normalized_data <- get_normalized_data(facs_raw_data, facs_meta_data,
                                       leukocytes, low_percent_count = 0.05,
                                       high_percent_count = 1,least_std = 1.2)
saveRDS(object=as.matrix(normalized_data$normalizded_data),
        file="leukocytes_normalized_data.rds")
saveRDS(object = normalized_data$meta_data, file = "leukocytes_meta_data.rds")
saveRDS(object = normalized_data$count_data, file = "leukocytes_count_data.rds")
