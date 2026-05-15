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

models <- sapply(norms, function(norm) { # Loading the results
  readRDS(paste0("scz_expression_rstan_regression_", norm, ".rds"))
})

pred <- setNames(
  lapply(norms, function(norm){
    ggpredict( # Using ggpredict to predict values from the model
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

invisible(
  lapply(norms, function(norm){ # Plotting in loop
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
      theme_minimal() +
      coord_cartesian(ylim = c(0, 15))
    print(plot)
  })
)
