library(doParallel)
library(optparse)
library(matrixStats)
library(caTools)
library(tidyverse)
library(here)
library(SIMPLEs)

## * load option
option_list <- list(
  make_option(c("--cnm"),
    action = "store",
    type = "integer", default = 1000
  )
)
args <- option_list %>%
  OptionParser(option_list = .) %>%
  parse_args()


load_sampled_data <- function(mycnm) {
  fnm <- paste0("leu_simple_cell_", mycnm, ".rds")
  readRDS(here::here("TabulaMuris", "supplement", "timer", fnm))
}


run_simple <- function() {
  logtpm_gbc <- load_sampled_data(args$cnm)
  SIMPLE(
    dat = logtpm_gbc, K0 = 10,
    M0 = 1, p_min = 0.4, cutoff = 0.1, min_gene = 1000,
    max_lambda = TRUE
  )
  return(0)
}

mytimer <- system.time({
  run_simple()
})
print(mytimer)

saveRDS(mytimer, file = here::here(
  "TabulaMuris", "supplement", "timer",
  paste0("simple_timer_", args$cnm, ".rds")
))
