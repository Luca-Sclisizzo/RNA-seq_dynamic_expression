suppressPackageStartupMessages({
library('lme4')
library('ggplot2')
library('dplyr')
library('splines')
library('ggeffects')
library('splineplot')
})

sample_metadata <- read.csv(
  "mRNA-seq_Sample metadata.csv",
  sep = ';',
  na.strings = c("NA", "")) %>%
  filter(!is.na(Window)) %>%
  filter(Ethnicity == 'European')
norms <- c('TMM', 'RLE', 'upperquartile')
models <- sapply(norms, function(norm) {
  readRDS(paste0("scz_expression_rstan_regression_", norm, ".rds"))
})

base <- model.frame(models[[1]])
newdata <- base[rep(1, 200), , drop = FALSE]
newdata$Days <- seq(
  min(sample_metadata$Days),
  max(sample_metadata$Days),
  length.out = 200
)

newdata$Sex <- factor("M", levels = levels(sample_metadata$Sex))
newdata$Sequencing.Site <- factor("YALE", levels = levels(sample_metadata$Sequencing.Site))

pred <- ggpredict(
  models[[3]],
  terms = "Days [0:13000 by=100]",
  condition = c(
    Sex = "M",
    Sequencing.Site = "YALE"
  )
)

plot(pred)

pred <- lapply(norms, function(norm){
  ggpredict(
    models[[norm]],
    terms = "Days [all]",
    condition = c(
      Sex = "M",
      Sequencing.Site = "YALE"
    )
  )
  
})
                  