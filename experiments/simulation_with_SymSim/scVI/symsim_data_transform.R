library(here)
library(tidyverse)

p <- list(nrep = 20, ncell = 300, ngene = 1000)

symsim_to_csv <- function(myseed) {
  to_csv <- function(tag = "umi") {
    if (tag == "umi") {
      obs_count <- obs_umi$counts
    } else if (tag == "nonumi") {
      obs_count <- obs_nonumi$counts
    } else {
      stop("only umi or nonumi.")
    }
    rownames(obs_count) <- seq_len(nrow(obs_count))
    colnames(obs_count) <- seq_len(ncol(obs_count))
    fnm <- stringr::str_c("sim", p$ncell, p$ngene, myseed, tag, ".csv", sep = "_")
    write.csv(obs_count, file = here::here("scVI", "data", "symsim", fnm))
  }
  fnm <- stringr::str_c("sim", p$ncell, p$ngene, myseed, ".RData", sep = "_")
  load(here::here("simutool", "jobs", "symsim_data", fnm))
  to_csv("umi")
  to_csv("nonumi")
}

seq_len(p$nrep) %>% map(symsim_to_csv)
