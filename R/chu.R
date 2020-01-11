#' Single cell RNASeq and bulk RNASeq data of human embryonic stem cell differentiation.
#'
#' Data from a study of human embryonic stem cell differentiation towards definitive endoderm. It includes both the single cell RNASeq and the bulk RNASeq data. Here we randomly choose part of the data.
#'
#' @docType data
#'
#' @usage data(chu)
#'
#' @format a list includes the following fields
#' \describe{
#'   \item{chu_normalized_data}{single cell RNASeq data, i.e., 
#'   3000 * 763 matrix with each row as gene and each column as cell.
#'   Each value is the log normalized gene expression.} 
#'   \item{chu_bulk_mean}{the average expression of the 3000 genes in the bulk RNASeq data} 
#'   \item{differential_genes}{3000 * 7 binary matrix. Each row is a gene, and each column is a unique cell type.}
#'   \item{chu_cell_type}{a character array with the length of 763.}
#'  }
#'
#' @keywords datasets
#'
#' @references Chu, L.F. et al., (2016) Genome Biology 17:1-20
#' (\href{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1033-x}{Genome Biology})
#'
#' @source \href{https://hemberg-lab.github.io/scRNA.seq.datasets/}{scRNA-Seq Datasets}
#'
#' @examples
#' data(chu)
#' single_cell_data <- chu$chu_normalized_data
#' bulk_rna_data  <- chu$chu_bulk_mean
#' cell_type <- chu$chu_cell_type
"chu"
