suppressPackageStartupMessages({
library('lme4')
library('ggplot2')
library('dplyr')
library('splines')
library('ggeffects')
library('splineplot')
library('bayesplot')
library('loo')
library('tidyr')
library('rstanarm')
})
norms <- c('TMM','RLE','upperquartile', 'none')


##### Frequentist fit ##### 
models <- sapply(norms, function(norm) { # Loading the results
  readRDS(paste0("scz_expression_freq_regression_", norm, ".rds"))
})


# Preparing a newdata df to predict on
newdata <- expand.grid(
  Window = seq(1, 9, length.out = 200),
  Sex = factor("F",levels = levels(model.frame(models$TMM)$Sex)),
  Sequencing.Site = factor("YALE",levels = levels(model.frame(models$TMM)$Sequencing.Site)),
  area = factor("DFC",levels = levels(model.frame(models$TMM)$area)))
# Predict
for(norm in norms){
  newdata[[paste0("pred_", norm)]] <- predict(
    models[[norm]],
    newdata = newdata,
    type = "response",
    re.form = NA
  )
} 
newdata <- pivot_longer(
  newdata,
  cols = starts_with("pred_"),
  names_to = "normalization",
  values_to = "prediction"
)
newdata$normalization <- sub("pred_", "", newdata$normalization)

# Plotting
windownames <- c("8-9pcw","12-13pcw","16-17pcw","19-22pcw","35pcw \n 4mos","0.5-2.5y","3-11y","13-19y","21-40y")

trajecotry_plot_frequentist <- ggplot(
  newdata,
  aes(x = Window, y = prediction, color = normalization, group = normalization)
) +
  geom_line(linewidth = 0.6) +
  scale_x_continuous(
    breaks = seq(1:9),
    labels = windownames
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  geom_vline(xintercept = 5, linetype = "dashed", col = "grey") +
  labs(
    x = "Developmental window",
    y = "Predicted expression"
  )
print(trajecotry_plot_frequentist)

fe <- fixef(models$TMM)
ci <- confint(models$TMM, method = "Wald")
ci_fe <- ci[names(fe), ]
df_coef <- data.frame(
  term = names(fe),
  estimate = as.numeric(fe),
  lower = ci_fe[,1],
  upper = ci_fe[,2]
) %>%
  dplyr::filter(!grepl("ns\\(Window", term))


ggplot(df_coef,
       aes(x = reorder(term, estimate),
           y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  coord_flip() +
  theme_minimal() +
  labs(x = "", y = "Effect size (β)")




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
