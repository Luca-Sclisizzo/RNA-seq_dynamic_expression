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
newdata$Window <- seq(
  min(sample_metadata$Window),
  max(sample_metadata$Window),
  length.out = 200
)

newdata$Sex <- factor("M", levels = levels(sample_metadata$Sex))
newdata$Sequencing.Site <- factor("YALE", levels = levels(sample_metadata$Sequencing.Site))

pred <- ggpredict(
  models[[2]],
  terms = "Window [all]",
  condition = c(
    Sex = "M",
    Sequencing.Site = "YALE"
  )
)
plot(pred)

pred <- setNames(
  lapply(norms, function(norm){
    ggpredict(
      models[[norm]],
      terms = "Window [all]",
      condition = c(
        Sex = "M",
        Sequencing.Site = "YALE"
      )
    )
  }),
  norms
)

lapply(norms, function(norm){
  df <- as.data.frame(pred[[norm]])
  plot <- ggplot(df, aes(x = x, y = predicted)) +
          geom_line(linewidth = 1) +
          geom_ribbon(
            aes(ymin = conf.low, ymax = conf.high),
            alpha = 0.2) +
          labs(
            title = paste("Norm:",norm),
            x = "Window",
            y = "Predicted value") +
          theme_minimal()
  print(plot)
})
