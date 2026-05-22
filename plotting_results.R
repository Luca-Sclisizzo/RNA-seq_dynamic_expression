suppressPackageStartupMessages({
library('lme4')
library('ggplot2')
library('dplyr')
library('splines')
library('ggeffects')
library('splineplot')
library('bayesplot')
library('loo')
library('rstanarm')
})
norms <- c('TMM','RLE','upperquartile')


##### Frequentist fit ##### 
models <- sapply(norms, function(norm) { # Loading the results
  readRDS(paste0("scz_expression_freq_regression_", norm, ".rds"))
})

newdata <- expand.grid(
  Window = seq(1,9,length.out = 200),
  Sex = "F",
  Sequencing.Site = 'YALE',
  area = "PFC"
)

pred <- setNames(
  lapply(norms, function(norm){
    newdata$pred <- predict(models[[norm]],
                            newdata = newdata,
                            type = "response",
                            re.form = NA)
}), norm
)




windownames <- c("8-9pcw","12-13pcw","16-17pcw","19-22pcw",
                 "35pcw \n 4mos","0.5-2.5y","3-11y","13-19y","21-40y")

invisible(
  lapply(norms, function(norm){ # Plotting in loop
    df <- as.data.frame(pred[[norm]])
    ggplot(df, aes(x = Window, y = pred, group = 1)) +
      geom_line() +
      scale_x_continuous(
        breaks = seq_along(windownames),
        labels = windownames
      )
  })
)




# pred <- setNames(
#   lapply(norms, function(norm){
#     ggpredict( # Using ggpredict to predict values from the model
#       models[[norm]],
#       bias_correction = TRUE, # suggested by the package
#       terms = "Window [all]",
#       condition = c(
#         Sex = "M",
#         Sequencing.Site = "YALE"
#       )
#     )
#   }),
#   norms
# )
# 
# invisible(
#   lapply(norms, function(norm){ # Plotting in loop
#     df <- as.data.frame(pred[[norm]])
#     plot <- ggplot(df, aes(x = x, y = predicted)) +
#       geom_line(linewidth = 1) +
#       geom_ribbon(
#         aes(ymin = conf.low, ymax = conf.high),
#         alpha = 0.2) +
#       labs(
#         title = paste("Norm:",norm),
#         x = "Window",
#         y = "Predicted value") +
#       theme_minimal() +
#       coord_cartesian(ylim = c(0, 15))
#     print(plot)
#   })
# )

##### Bayesian fit ##### 
TMM_meanfied_model <- readRDS('scz_expression_bayes_regression_meanfield_TMM.rds')
posterior <- as.matrix(TMM_meanfied_model)
posterior <- posterior[, startsWith(colnames(posterior), "ns")]
colnames(posterior) <- c('1','2','3','4')

posterior <- as.data.frame(TMM_meanfied_model)
posterior <- draws[, startsWith(colnames(draws), "ns")]
mcmc_intervals(posterior, prob = 0.80)



colnames(ns_cols) <- c('1','2','3','4')




mcmc_intervals(ns_cols, prob = 0.80)
mcmc_areas(posterior,
               prob = 0.50, prob_outer = 0.95,
           pars = c('ns(Days, df = 4)1', 'ns(Days, df = 4)2', 'ns(Days, df = 4)3', 'ns(Days, df = 4)4'))


loo_TMM <- loo(TMM_meanfied_model, save_psis = TRUE)
plot(loo_TMM)
rstan::get_stanmodel(TMM_meanfied_model$stanfit)
