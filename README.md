# SIMPLEs: single-cell RNA sequencing imputation and cell clustering methods by modeling gene module variation

Authors:  Zhirui Hu, Songpeng Zu

## Introduction
SIMPLE is a R package that imputes "dropout" values (i.e. zero or small entries)
in the single cell RNASeq data based on cell similarities and gene correlations.
The imputed matrix can be used for dimension reduction and visualization of the
cell spectrum, to identify markers between different samples and construct gene
co-expression network. 

### How it works
SIMPLEs iteratively clusters the cells, identifies
correlated gene modules, infers the probability of the dropout event for each
zero entry, and imputes zeros within each cluster utilizing the expression of other
correlated genes. It will impute the technical zeros (or "dropout") while
maintain the biological zeros at low expression level.

The imputation process is based on the correlations between genes within similar
cell types, which is modeled by several common latent factors or gene modules as
well as the gene-specific dropout rate. Although the dropout rate can be
estimated from the empirical distribution of gene expression in the scRNASeq, it
could interference with estimating the gene correlation structure, especially
for lowly expressed genes. 

Integrating with the corresponding bulk RNASeq data
which serves as the average gene expression across cells, provides extra sources
of information on the dropout rate per gene. It can give a better estimate of
the "dropout" rate which would influence how much the data should be imputed. 

We called our method integrating bulk RNASeq data as **SIMPLE-B**, otherwise
**SIMPLE**. We referred our toolset including SIMPLE and SIMPLE-B as
**SIMPLEs**.

## Installation
The package is not on CRAN yet. You can use the following codes in `R` to
install it.

```r
#> use devtools package to install the repo from github
#> - if you want to follow the latest but unstable version, you can set 
#>   the ref as"develop".
#> - if you want to build the vignettes for reference, you can set 
#>   build_vignettes as TRUE.
devtools::install_github("xyz111131/SIMPLEs", ref="master", build_vignettes=FALSE)
```
## Vignettes
You can view the vignettes in the package as reference.
```r
vignettes(topic="SIMPLEs_example", package="SIMPLEs")
```

## Data
In the vignettes, we introduced how to simulate the single cell RNASeq data.
Furthermore, we provided a real dataset from a recent study. You can view
the details of the data, and run the analyses using SIMPLEs.
```r
# describe the details of the real dataset
?SIMPLEs::chu
```

