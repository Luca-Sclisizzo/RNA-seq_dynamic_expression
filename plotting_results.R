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
# Trajectory plots
newdata_bayes <- expand.grid(
  Window = seq(1, 9, length.out = 200),
  Sex = factor("F",levels = levels(model.frame(model_single_gene_first_score)$Sex)),
  Sequencing.Site = factor("YALE",levels = levels(model.frame(model_single_gene_first_score)$Sequencing.Site)),
  area = factor("DFC",levels = levels(model.frame(model_single_gene_first_score)$area))
)
pred_mat_1gene <- rstanarm::posterior_linpred(
  model_single_gene_first_score,
  newdata = newdata_bayes,
  re.form = NA
)
df_bayes <- data.frame(
  Window = newdata_bayes$Window,
  pred_first_gene_mean = colMeans(pred_mat_1gene),
  pred_first_gene_low  = apply(pred_mat_1gene, 2, quantile, probs = 0.025),
  pred_first_gene_high = apply(pred_mat_1gene, 2, quantile, probs = 0.975)
)


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
#gene_50_spline_areas <- bayesplot::mcmc_areas(model_single_gene_50_score, pars = par)

# Trajectory plots
pred_mat_50Gene <- rstanarm::posterior_linpred(
  model_single_gene_50_score,
  newdata = newdata_bayes,
  re.form = NA
)
df_bayes <- df_bayes %>%
  mutate(
    pred_50_gene_mean = colMeans(pred_mat_50Gene),
    pred_50_gene_low  = apply(pred_mat_50Gene, 2, quantile, probs = 0.025),
    pred_50_gene_high = apply(pred_mat_50Gene, 2, quantile, probs = 0.975)
    )
df_bayes_long <- df_bayes %>%
  pivot_longer(
    cols = -Window,
    names_to = c("gene", "stat"),
    names_pattern = "pred_(.*)_(mean|low|high)",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  )
# Plotting with 50% ribbon
spline_traj_single_gene <- ggplot(df_bayes_long, aes(x = Window, y = mean, color = gene, fill = gene)) +
                                  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, color = NA) +
                                  geom_line(linewidth = 0.8) +
                                  theme_classic(base_size = 13) +
                                  scale_x_continuous(
                                    breaks = 1:9,
                                    labels = windownames
                                  ) +
                                  theme_classic(base_size = 13) +
                                  theme(
                                    axis.text.x = element_text(angle = 45, hjust = 1),
                                    legend.title = element_blank()
                                  ) +
                                  geom_vline(xintercept = 5, linetype = "dashed", color = "grey") +
                                  labs(
                                    x = "Developmental window",
                                    y = "Predicted expression (posterior median)",
                                    color = "Gene",
                                    fill = "Gene",
                                    caption = "Light shaded areas represent 50% credible intervals (25–75% posterior quantiles)"
                                  )

# Clusterwise Models & plots ----------------------------------------------
#### Clusterwise fixed interaction (ns(Window)*Modules) ####
model_clusterwise_interaction <- readRDS("~/RNA-seq_training/scz_expression_bayes_regression_MCMC_rawcounts_clusterwise_interactions_TMM.rds")
yrep_clusterwise_interaction <- rstanarm::posterior_predict(model_clusterwise_interaction, draws = 500)
median_clusterwise_interaction <- ppc_stat(model_clusterwise_interaction$y, yrep_clusterwise_interaction, stat = "median", discrete = TRUE) + ggtitle('Median') + 
  theme(plot.title = element_text(hjust = 0.5, size = 12)) + theme(legend.position = "none")
sd_clusterwise_interaction <- ppc_stat(model_clusterwise_interaction$y, yrep_clusterwise_interaction, stat = "sd", discrete = TRUE) + ggtitle('Sd') + 
  theme(plot.title = element_text(hjust = 0.5, size = 12)) + theme(legend.position = "none")
sum_0s_clusterwise_interaction <- ppc_stat(model_clusterwise_interaction$y, yrep_clusterwise_interaction, stat=function(x) sum(x==0)) + xlim(0, 5) + ggtitle('Sum of 0s') +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

#### Clusterwise interaction  w gene variability(ns(Window)*Modules) ####
# model_clusterwise_interaction_w_variability <- readRDS('/Users/lucasclisizzo/RNA-seq_training/scz_expression_bayes_regression_rstan_clusterwise_gene_variability_TMM.rds')
# model <- cmdstan_model("rstan_model.stan") # model with y_rep generation
# draws_df <- model_clusterwise_interaction_w_variability$draws(
#    format = "draw.array"
#  )[1:100, , ]  # prime 100 draws
# mcmc_trace(draws_df, pars = c("mu", "sigma"))
# par <- c(
#   "alpha",
#   as.vector(outer(1:9, 1:4, \(i, j) paste0("beta_module[", i, ",", j, "]")))
# )
# model_clusterwise_interaction_w_variability$summary(variables = par) |>
#   dplyr::filter(rhat > 1.01 | ess_bulk < 400)

# gq_fit <- model$generate_quantities(
#   fitted_params = draws_df,
#   data = stan_data
# )
# y_rep <- gq_fit$draws("y_rep", format = "matrix")  # [S x N]
# 
# 
# ppc_stat(stan_data$y, y_rep, stat = "median", discrete = TRUE) + ggtitle('Median') + 
#   theme(plot.title = element_text(hjust = 0.5, size = 12)) + theme(legend.position = "none")




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
