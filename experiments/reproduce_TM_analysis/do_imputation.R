import::from(here, here)
library(doParallel)

## * load comparison methods
import::from(SIMPLEs, SIMPLE)
import::from(optparse, make_option, parse_args, OptionParser)

## * setup parser

set_parser <- function() {
  option_list  <- list(
    make_option("--M0",action="store", type="integer", default=1, help="number of init cluster number"),
    make_option(c("-p", "--pm"), action="store", type="numeric", default=0.3,
                help="minumum ratio of expressed counts for the high quality genes selection, we'd better keep the high quality genes number at least 2000 or 3000."),
    make_option(c("-c", "--cutoff"), action="store", type="numeric", default=0.1,
                help="values below this cutoff are treated as dropout."),
    make_option("--mingene", action="store", type="integer",default=1000,
                help="minimun number of genes of high quality. It's used when no more high quality are chosen than this."),
    make_option("--cores", action="store", type="integer",default=8,
                help="number of cores to run parallelly")
  )
  return(parse_args(OptionParser(option_list=option_list)))
}

## * run imputation
impute_with_SIMPLE  <- function(data, K0=10, M0=1, pm = 0.3, cutoff=0.1, mingene=1000, max_lambda=T){
  simple_result <- SIMPLE(dat=as.matrix(data), K0=K0, M0=M0, p_min = pm, cutoff=cutoff,min_gene = mingene, max_lambda=max_lambda)
  save_file_name <- paste0("SIMPLE_","K0-",K0, "_M0-", M0, "_pm-", pm,"_cufoff-",cutoff, "_maxlambda_result.rds")
  saveRDS(object=simple_result, file=save_file_name)
}

## * main
normalized_data <- readRDS(file=here::here("leukocytes_normalized_data.rds"))
parser <- set_parser()
registerDoParallel(cores = parser$cores)

message("IMPUTATION BEGINS with method ", parser$method)
impute_with_SIMPLE(data=normalized_data, K0=parser$K0, M0=parser$M0, pm=parser$pm,
                     cutoff=parser$cutoff, mingene=parser$mingene, max_lambda=T)
message("IMPUTATION END.")
