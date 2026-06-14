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
library('patchwork')
})
norms <- c('TMM','RLE','upperquartile', 'none')
windownames <- c("8-9pcw","12-13pcw","16-17pcw","19-22pcw","35pcw \n 4mos","0.5-2.5y","3-11y","13-19y","21-40y")
par <- c('ns(Window, df = 4)1', 'ns(Window, df = 4)2', 'ns(Window, df = 4)3', 'ns(Window, df = 4)4')

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
#print(trajecotry_plot_frequentist)
# 
# fe <- fixef(models_freq$TMM)
# ci <- confint(models_freq$TMM, method = "Wald")
# ci_fe <- ci[names(fe), ]
# df_coef <- data.frame(
#   term = names(fe),
#   estimate = as.numeric(fe),
#   lower = ci_fe[,1],
#   upper = ci_fe[,2]
# ) %>%
#   dplyr::filter(!grepl("ns\\(Window", term))
# 
# 
# ggplot(df_coef,
#        aes(x = reorder(term, estimate),
#            y = estimate)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
#   coord_flip() +
#   theme_minimal() +
#   labs(x = "", y = "Effect size (β)")


# Bayesian Complete pooling & plots ----------------------------------------------------
#### Model with pseudocounts ####
models_bayes <- sapply(norms, function(norm) { # Loading the results
  readRDS(paste0("scz_expression_bayes_regression_MCMC_", norm, ".rds"))
}, simplify = FALSE)

posteriors_bayes <- lapply(norms, function(norm) {
  as.array(models_bayes[[norm]])
})
names(posteriors_bayes) <- norms
for (norm in norms) { # Renaming the spline coeffiecients names
  dimnames(posteriors_bayes[[norm]])[[3]][2:5] <- paste0("spline_df", 1:4)
}
### Creating posterior plots using Bayesplot
mcmc_areas_plots <- list()
mcmc_dens_overlay_plots <- list()

for (norm in norms){
  pars <- c('spline_df1','spline_df2','spline_df3','spline_df4')
  
  mcmc_areas_plots[[norm]] <- mcmc_areas(posteriors_bayes[[norm]], 
                                         par = pars) + ggtitle(norm) + 
                                        theme(plot.title = element_text(hjust = 0.5, size = 10))
  mcmc_dens_overlay_plots[[norm]] <- mcmc_dens_overlay(posteriors_bayes[[norm]],
                                                       par = pars) + ggtitle(norm) + 
                                                      theme(plot.title = element_text(hjust = 0.5, size = 10))
  if (norm %in% c('RLE', 'none')){ # Removing y lables and ticks for the wrap
    mcmc_areas_plots[[norm]] <- mcmc_areas_plots[[norm]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
}

### Posterior predictive checks
yrep_bayes <- setNames(
  lapply(norms, function(norm){
    yrep_bayes <- rstanarm::posterior_predict(models_bayes[[norm]], draws = 500)
  }), norms)
# ppc_plots <- lapply(norms, function(norm) {
#   ppc_stat(
#     y = models_bayes[[norm]]$y,
#     yrep = yrep_bayes[[norm]],
#     stat = "median",
#     discrete = TRUE
#   )
# })
# names(ppc_plots) <- norms

### Trajectory plotting
newdata_bayes <- expand.grid( # Creating a dataset of new values
  Window = seq(1, 9, length.out = 200),
  Sex = 'F',
  Sequencing.Site = 'YALE',
  area = 'DFC'
  )

Y_hat_list <- list()
for (norm in norms){ # Using the posterior predictive on the 'holdout' data  
  Y_hat <- posterior_epred(
    models_bayes[[norm]],
    newdata = newdata_bayes,
    re.form = NA
  )
  Y_hat_list[[norm]] <- Y_hat
}
# Creating posterior predictive point estimates
# I will summarize the posterior predictive using the median and 50% Credibility Interval (quantile 25% - 75%)
summary_Y_pred_list <- lapply(Y_hat_list, function(Y_hat){
  data.frame(
    Y_med = apply(Y_hat, 2, median), # Median
    Y_lb  = apply(Y_hat, 2, quantile, probs = 0.25), # Lower Bound (25%)
    Y_ub  = apply(Y_hat, 2, quantile, probs = 0.75) # Upper Bound (75%)
  )
})
for (norm in names(summary_Y_pred_list)) { # Adding the Posterior Predictive summaries to the newdata_bayes df
  tmp <- summary_Y_pred_list[[norm]]
  newdata_bayes[[paste0("Y_med_", norm)]] <- tmp$Y_med
  newdata_bayes[[paste0("Y_lb_",  norm)]] <- tmp$Y_lb
  newdata_bayes[[paste0("Y_ub_",  norm)]] <- tmp$Y_ub
}
newdata_bayes <- newdata_bayes %>% # Transforming the newdata_bayes to long in order to plot it
  pivot_longer(
    cols = matches("^Y_(med|lb|ub)_"),
    names_to = c("stat", "norm"),
    names_pattern = "Y_(med|lb|ub)_(.*)",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  )

# Plotting with ribbon (50% CI)
trajectory_plot_bayesian <- ggplot(newdata_bayes, aes(x = Window, y = med, color = norm, fill = norm)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.2, color = NA) + # Credibility interval 50%
  geom_line(linewidth = 0.8) +
  scale_x_continuous(
    breaks = 1:9,
    labels = windownames
  ) +
  theme_classic(base_size = 13) +
  facet_wrap(~norm) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "grey") +
  labs(
    x = "Developmental window",
    y = "Predicted expression (posterior median)",
    caption = "Note: shaded regions represent 50% credible intervals (25–75% posterior quantiles)"
  )
#print(trajectory_plot_bayesian)

#### Model with RawCounts ####
model_bayes_raw_counts_TMM <- readRDS('/Users/lucasclisizzo/RNA-seq_training/scz_expression_bayes_regression_MCMC_rawcounts_TMM.rds')
yrep_bayes_raw_counts_TMM <- rstanarm::posterior_predict(model_bayes_raw_counts_TMM, draws = 500)
median_PPC <-   ppc_stat(y = model_bayes_raw_counts_TMM$y, yrep = yrep_bayes_raw_counts_TMM, stat = "median",discrete = TRUE) + ggtitle('Median') +
  theme(plot.title = element_text(hjust = 0.5, size = 12)) + theme(legend.position = "none")
sd_PPC <-   ppc_stat(y = model_bayes_raw_counts_TMM$y, yrep = yrep_bayes_raw_counts_TMM, stat = "sd",discrete = TRUE) + ggtitle('Sd') +
  theme(plot.title = element_text(hjust = 0.5, size = 12)) + theme(legend.position = "none")
sum_0s <- ppc_stat(model_bayes_raw_counts_TMM$y, yrep_bayes_raw_counts_TMM, stat=function(x) sum(x==0)) + ggtitle('Sum of 0s') + 
  theme(plot.title = element_text(hjust = 0.5, size = 12))

# Single Gene Models & plots ------------------------------------------------------------
####  Single gene model: first gene as score ####
model_single_gene_first_score <- readRDS('/Users/lucasclisizzo/RNA-seq_training/scz_expression_bayes_regression_MCMC_rawcounts_single_gene_TMM.rds')
yrep_model_gene_first <- rstanarm::posterior_predict(model_single_gene_first_score, draws = 500)
median_single_first_gene <- ppc_stat(model_single_gene_first_score$y, yrep_model_gene_first, stat = "median", discrete = TRUE) + ggtitle('Median') + 
  theme(plot.title = element_text(hjust = 0.5, size = 12)) + theme(legend.position = "none")
sd_single_first_gene <- ppc_stat(model_single_gene_first_score$y, yrep_model_gene_first, stat = "sd", discrete = TRUE) + ggtitle('Sd') + 
  theme(plot.title = element_text(hjust = 0.5, size = 12)) + theme(legend.position = "none")
sum_0s_single_first_gene <- ppc_stat(model_single_gene_first_score$y, yrep_model_gene_first, stat=function(x) sum(x==0)) + xlim(0, 5) + ggtitle('Sum of 0s') +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

# Spline coefficients distributions
gene_first_spline_areas <- bayesplot::mcmc_areas(model_single_gene_first_score, pars = par)

#### Single gene model: 50th gene as score ####
model_single_gene_50_score <- readRDS('/Users/lucasclisizzo/RNA-seq_training/scz_expression_bayes_regression_MCMC_rawcounts_single_gene_50_th_gene_TMM.rds')
yrep_model_gene_50 <- rstanarm::posterior_predict(model_single_gene_50_score, draws = 500)
median_single_50_gene <- ppc_stat(model_single_gene_50_score$y, yrep_model_gene_50, stat = "median", discrete = TRUE) + ggtitle('Median') + 
  theme(plot.title = element_text(hjust = 0.5, size = 12)) + theme(legend.position = "none")
sd_single_50_gene <- ppc_stat(model_single_gene_50_score$y, yrep_model_gene_50, stat = "sd", discrete = TRUE) + ggtitle('Sd') + 
  theme(plot.title = element_text(hjust = 0.5, size = 12)) + theme(legend.position = "none")
sum_0s_single_50_gene <- ppc_stat(model_single_gene_50_score$y, yrep_model_gene_50, stat=function(x) sum(x==0)) + xlim(0, 5) + ggtitle('Sum of 0s') +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

# Spline coefficients distributions
gene_50_spline_areas <- bayesplot::mcmc_areas(model_single_gene_50_score, pars = par)

patchwork::wrap_plots(list(gene_first_spline_areas, gene_50_spline_areas))
# Clusterwise Models & plots ----------------------------------------------
####  Clusterwise random intercept #### 
model_clusterwise_intercept <- readRDS("~/RNA-seq_training/scz_expression_bayes_regression_MCMC_rawcounts_clusterwise_TMM.rds")
yrep_clusterwise_intercept <- rstanarm::posterior_predict(model_clusterwise_intercept, draws = 500)
median_clusterwise_intercept <- ppc_stat(model_clusterwise_intercept$y, yrep_clusterwise_intercept, stat = "median", discrete = TRUE) + ggtitle('Median') + 
  theme(plot.title = element_text(hjust = 0.5, size = 12)) + theme(legend.position = "none")
sd__clusterwise_intercept <- ppc_stat(model_clusterwise_intercept$y, yrep_clusterwise_intercept, stat = "sd", discrete = TRUE) + ggtitle('Sd') + 
  theme(plot.title = element_text(hjust = 0.5, size = 12)) + theme(legend.position = "none")
sum_0s__clusterwise_intercept <- ppc_stat(model_clusterwise_intercept$y, yrep_clusterwise_intercept, stat=function(x) sum(x==0)) + xlim(0, 5) + ggtitle('Sum of 0s') +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

#### Clusterwise fixed interaction (ns(Window)*Modules) ####
model_clusterwise_interaction <- readRDS("~/RNA-seq_training/scz_expression_bayes_regression_MCMC_rawcounts_clusterwise_interactions_TMM.rds")
yrep_clusterwise_interaction <- rstanarm::posterior_predict(model_clusterwise_interaction, draws = 500)
median_single_50_gene <- ppc_stat(model_clusterwise_interaction$y, yrep_clusterwise_interaction, stat = "median", discrete = TRUE) + ggtitle('Median') + 
  theme(plot.title = element_text(hjust = 0.5, size = 12)) + theme(legend.position = "none")
sd_single_50_gene <- ppc_stat(model_clusterwise_interaction$y, yrep_clusterwise_interaction, stat = "sd", discrete = TRUE) + ggtitle('Sd') + 
  theme(plot.title = element_text(hjust = 0.5, size = 12)) + theme(legend.position = "none")
sum_0s_single_50_gene <- ppc_stat(model_clusterwise_interaction$y, yrep_clusterwise_interaction, stat=function(x) sum(x==0)) + xlim(0, 5) + ggtitle('Sum of 0s') +
  theme(plot.title = element_text(hjust = 0.5, size = 12))



### Model comparison (eventually)
# loo_activate <- FALSE
# if (!loo_activate) {
#   warning(
#     "LOO is disabled: Leave-One-Out will not be computed due to computational cost.\n",
#     "Set loo_activate <- TRUE to enable it."
#   )
# } else{
#   cores <- 4 # Specified in case loo() needs to partial refit the model removing the influential observations
#   loo_bayes <- lapply(norms, function(norm) { # Creating posterior infomation criterias
#     loo(models_bayes[[norm]], k_threshold = 0.7)
#   })
#   names(loo_bayes) <- paste0("loo_", norms)
# }
# loo_compare <- loo::loo_compare(loo_bayes[["loo_TMM"]], loo_bayes[["loo_RLE"]], loo_bayes[["loo_upperquartile"]], loo_bayes[["loo_none"]])
