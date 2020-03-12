# Reproduction of the analysis on the Tabula Muris database.

Author: Songpeng Zu

## 1. Prepare the packages and dataset.
* Run the script *prepare_packages.R* to install the R packages needed in this
analysis.

* Download Tabula Muris dataset
You can download dataset directly from the original [Tabula Muris
repository](https://github.com/czbiohub/tabula-muris). Here we write a Makefile
to simplify the download process. our analysis focuses on FACS-based full length
transcript analysis data.

  * Download TM_fact_mat.rds (518M) and TM_facs_metadata.csv (696K) to the current directory.
    ```shell
    make fact_all
    ```
   * According to the [Tabula Muris vignettes
     repository](https://github.com/czbiohub/tabula-muris-vignettes), we use
     [Cell Ontology web](http://obofoundry.org/ontology/cl.html) to select the
     immune cells. You can download cell_ontology.obo (888K) to the current
     directory.
     ```shell
     make cell_ontology
     ```

## 2. Data preprocessing
* Run the script *data_preprocessing.R* to do data processing including
  normalization, select leukocyte cells, and so on. This script will generate
  three data named *leukocytes_normalizded_data.rds*,
  *leukocytes_meta_data.rds*, *leukocytes_count_data.rds*.

## 3. Do imputation
Use SIMPLEs to impute the leukocytes single cell data above. The meanings of
different parameters can be found by using `?SIMPLE` after loading our packages.
`library(SIMPLEs)`.
```
Rscript --vanilla do_imputation.R -k 10 --M0 1 -p 0.3 -c 0.1 --cores 8 
```
This process might take several hours. In our case, it needs about 6-8 hours on
a 8-core single node. This will generate the file named
*SIMPLE_K0-10_M0-1_pm-0.3_cufoff-0.1_maxlambda_result.rds*. The file name will
be changed according to different parameters.

In our paper, Figure 8 shows the results of 20 different combinations of K (5,
10, 15, 20) and M (1,2,8,12,22). Note that pm is fixed as 0.4.

## 4. TSNE results for different parameters of SIMPLEs. 
* We use the function `handle_raw_data` in the script `analyze_imputation.R` to
get the Figure 8(b), Figure S11 and Figure S12.

* Then we use the function `handle_simple_tsne` and `handle_simple_draw_tsne`
  in the same script to get Figure 8(a).

## 5. Boxplot and dotplot for gene module analysis given M as 1 and K as 10.
* We use function `draw_simple_boxplot` in the script `analyze_imputation.R` to
  get the Figure 9(a) and Figure S13.

* Start from Line 387 till the end in the same script, we draw the dotplot for
  Figure 9(b). 
