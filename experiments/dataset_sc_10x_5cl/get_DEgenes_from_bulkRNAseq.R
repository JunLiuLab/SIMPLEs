library(DESeq2)
library(org.Hs.eg.db)
library(here)
library(tidyverse)
## * meta data
## infomation is summarized in the file GSE86337_family.soft or
## GSE86337_series_matrix.txt

## GSE86337_reverse.stranded.unfiltered.count.matrix.txt the matrix is tab
# seprated.
# row labels # "ID_REF" "GSM2300446" "GSM2300447" "GSM2300448"
# "GSM2300449" "GSM2300450" "GSM2300451" "GSM2300452" "GSM2300453" "GSM2300454"
# "GSM2300455"
# Sample_geo_accession "GSM2300446" "GSM2300447" "GSM2300448"
# "GSM2300449" "GSM2300450" "GSM2300451" "GSM2300452" "GSM2300453" "GSM2300454"
# "GSM2300455"
# Sample_title "A549 rep1" "A549 rep2" "H1975 rep1" "H1975 rep2"
# "HCC827 rep1" "HCC827 rep2" "H2228 rep1" "H2228 rep2" "H838 rep1" "H838 rep2"

## Gene Entrez is used here.
## Genome_build: hg19 for sample data processing

## * load data
orig_bulk_rnas <- read.csv(
  file = here::here(
    "10xGenomics", "bulk_RNAseq",
    "GSE86337_reverse.stranded.unfiltered.count.matrix.txt"
  ),
  header = TRUE, sep = "\t", row.names = 1
)

## * use gene symbol to align with scRNAseq.

orig_genesymbols <- AnnotationDbi::select(org.Hs.eg.db,
                                          keys = rownames(orig_bulk_rnas),
                                          column = "SYMBOL", keytype = "ENTREZID"
                                          )
## remove unmapped gene symbole rows
kept_rows <- which(!is.na(orig_genesymbols$SYMBOL))
bulk_rnas <- orig_bulk_rnas[kept_rows, ]
## remapping
genesymbols <- orig_genesymbols[kept_rows, ]$SYMBOL
rownames(bulk_rnas) <- genesymbols

colnames(bulk_rnas) <- c(
  rep("A549", 2), rep("H1975", 2), rep("HCC827", 2),
  rep("H2228", 2), rep("H838", 2)
)

## * DEseq2 analysis
groups <- factor(colnames(bulk_rnas))
dds <- DESeqDataSetFromMatrix(
  countData = bulk_rnas, colData = DataFrame(groups),
  design = ~groups
) %>% DESeq()
# use resultsNames(degs) to view the result names
## h1va means h1975 vs a549
dh1va <- results(dds, name = "groups_H1975_vs_A549")
h1va <- data.frame(
  log2fc = dh1va$log2FoldChange,
  padj = dh1va$padj,
  row.names = rownames(bulk_rnas)
)
h1va$padj[is.na(h1va$padj)] <- 1

dh2va <- results(dds, name = "groups_H2228_vs_A549")
h2va <- data.frame(
  log2fc = dh2va$log2FoldChange,
  padj = dh2va$padj,
  row.names = rownames(bulk_rnas)
)
h2va$padj[is.na(h2va$padj)] <- 1

dh8va <- results(dds, name = "groups_H838_vs_A549")
h8va <- data.frame(
  log2fc = dh8va$log2FoldChange,
  padj = dh8va$padj,
  row.names = rownames(bulk_rnas)
)
h8va$padj[is.na(h8va$padj)] <- 1


dhcva <- results(dds, name = "groups_HCC827_vs_A549")
hcva <- data.frame(
  log2fc = dhcva$log2FoldChange,
  padj = dhcva$padj,
  row.names = rownames(bulk_rnas)
)
hcva$padj[is.na(hcva$padj)] <- 1

## set degene demo: h1va[(abs(deseq_h1va$log2FoldChange) < 0.5)
##      & (deseq_h1va$padj > 0.05)]  <-  0

## * save results
fres <- list(
  gmap = genesymbols,
  h1975va549 = h1va,
  h2228va549 = h2va,
  h838va549 = h8va,
  hcc827va549 = hcva
)

saveRDS(object = fres, file = here::here("10xGenomics", "bulk_RNAseq", "DESeq2_result.RDS"))
