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

# Frequentist fit & plots -------------------------------------------------
models_freq <- sapply(norms, function(norm) { # Loading the results
  readRDS(paste0("scz_expression_freq_regression_", norm, ".rds"))
}, simplify = FALSE)


# Preparing a newdata df to predict on
newdata <- expand.grid(
  Window = seq(1, 9, length.out = 200),
  Sex = factor("F",levels = levels(model.frame(models_freq$TMM)$Sex)),
  Sequencing.Site = factor("YALE",levels = levels(model.frame(models_freq$TMM)$Sequencing.Site)),
  area = factor("DFC",levels = levels(model.frame(models_freq$TMM)$area)))
# Predict
for(norm in norms){
  newdata[[paste0("pred_", norm)]] <- predict(
    models_freq[[norm]],
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

fe <- fixef(models_freq$TMM)
ci <- confint(models_freq$TMM, method = "Wald")
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





# Bayesian fit & plots ----------------------------------------------------
models_bayes <- sapply(norms, function(norm) { # Loading the results
  readRDS(paste0("scz_expression_bayes_regression_MCMC_", norm, ".rds"))
}, simplify = FALSE)

posteriors_bayes <- sapply(norms, function(norm) { # Creating posteriors
  assign(
    paste0('posterior_',norm),
    as.array(models_bayes[[norm]])
    )
}, simplify = FALSE)

for (norm in norms) {
  names(posteriors_bayes[[norm]])[2:5] <- paste0("spline_df", 1:4)
}
invisible(
  lapply(norms, function(norm){
    pars <- c('ns(Window, df = 4)1','ns(Window, df = 4)2','ns(Window, df = 4)3','ns(Window, df = 4)4')
    title_results <- ggtitle(paste0('Resuls_',norm))
    print(
      mcmc_areas(posteriors_bayes[[norm]], 
                 par = pars) + title_results
    )
    print(
      mcmc_dens_overlay(posteriors_bayes[[norm]], 
                        par = pars) + title_results
    )
  })
)

loo_activate <- FALSE
if (!loo_activate) {
  warning(
    "LOO is disabled: Leave-One-Out will not be computed due to computational cost.\n",
    "Set loo_activate <- TRUE to enable it."
  )
} else{
  loo_bayes <- lapply(norms, function(norm) { # Creating posterior infomation criterias
    loo(models_bayes[[norm]])
  })
  names(loo_bayes) <- paste0("loo_", norms)
}




# Old code ----------------------------------------------------------------

# TMM_meanfied_model <- readRDS('scz_expression_bayes_regression_meanfield_TMM.rds')
# posterior <- as.data.frame(TMM_meanfied_model)
# posterior <- posterior[, startsWith(colnames(posterior), "ns")]
# colnames(posterior) <- c('1','2','3','4')
# 
# posterior <- as.data.frame(TMM_meanfied_model)
# posterior <- posterior[, startsWith(colnames(posterior), "ns")]
# mcmc_areas(posterior, prob = 0.80)
# 
# yrep <- posterior_predict(TMM_meanfied_model, draws = 500)
# ppc_dens_overlay(y = TMM_meanfied_model$y, 
#                  yrep = yrep)

# 
# color_scheme_set("brightblue")
# TMM_meanfied_model %>%
#   ppc_stat(y = TMM_meanfied_model$y,
#                    yrep = yrep,
#                    stat = "median")
# 
# count_zeros <- function(x) {sum(x == 0)}
# TMM_meanfied_model %>%
#   ppc_stat(y = TMM_meanfied_model$y,
#            yrep = yrep,
#            stat = count_zeros)
