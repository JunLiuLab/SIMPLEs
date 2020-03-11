# Reproduction of the analysis on the Tabula Muris database.

Authors: Songpeng Zu

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
  the data named *leukocytes_normalizded_data.rds*
