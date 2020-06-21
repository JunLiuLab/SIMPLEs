library(tidyverse)
library(here)


## * load original dataset.
logtpm_gbc <- readRDS(here::here(
  "TabulaMuris", "data",
  "leukocytes_normalized_data.rds"
))
cnm <- dim(logtpm_gbc)[2]

save_sampled_data <- function(mycnm) {
  fnm <- paste0("leu_simple_cell_", mycnm, ".rds")
  logtpm <- logtpm_gbc[, sample(1:cnm, mycnm)]
  saveRDS(object = logtpm, file = here::here("TabulaMuris", "supplement", "timer", fnm))
}

c(1000, seq(2000, 12000, 2000)) %>% map(save_sampled_data)
