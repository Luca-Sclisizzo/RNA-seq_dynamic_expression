suppressPackageStartupMessages({
library('lme4')
library('ggplot2')
library('dplyr')
library('splines')
library('ggeffects')
library('splineplot')
library('bayesplot')
library('loo')
})

##### Frequentist fit ##### 
models <- sapply(norms, function(norm) { # Loading the results
  readRDS(paste0("scz_expression_freq_regression_raw_counts_", norm, ".rds"))
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

##### Bayesian fit ##### 
TMM_meanfied_model <- readRDS('scz_expression_bayes_regression_meanfield_TMM.rds')
loo_TMM <- loo(TMM_meanfied_model, save_psis = TRUE)
plot(loo_TMM)
rstan::get_stanmodel(TMM_meanfied_model$stanfit)
