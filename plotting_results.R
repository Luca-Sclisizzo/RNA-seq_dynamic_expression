suppressPackageStartupMessages({
library('lme4')
library('ggplot2')
library('dplyr')
})

sample_metadata <- read.csv(
  "mRNA-seq_Sample metadata.csv",
  sep = ';',
  na.strings = c("NA", "")
)
norms <- c('TMM', 'RLE', 'upperquartile')
models <- sapply(norms, function(norm) {
  readRDS(paste0("scz_expression_rstan_regression_", norm, ".rds"))
})



