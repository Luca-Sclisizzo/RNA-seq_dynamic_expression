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

# Bayesian Complete pooling & plots ----------------------------------------------------
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
####  Single gene model: 1st gene as score ####
model_single_gene_first_score <- readRDS('/Users/lucasclisizzo/RNA-seq_training/scz_expression_bayes_regression_MCMC_rawcounts_single_gene_TMM.rds')
yrep_model_gene_first <- rstanarm::posterior_predict(model_single_gene_first_score, draws = 500)
median_single_first_gene <- ppc_stat(model_single_gene_first_score$y, yrep_model_gene_first, stat = "median", discrete = TRUE) + ggtitle('Median') + 
  theme(plot.title = element_text(hjust = 0.5, size = 12)) + theme(legend.position = "none")
sd_single_first_gene <- ppc_stat(model_single_gene_first_score$y, yrep_model_gene_first, stat = "sd", discrete = TRUE) + ggtitle('Sd') + 
  theme(plot.title = element_text(hjust = 0.5, size = 12)) + theme(legend.position = "none")
sum_0s_single_first_gene <- ppc_stat(model_single_gene_first_score$y, yrep_model_gene_first, stat=function(x) sum(x==0)) + ggtitle('Sum of 0s') +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

# Spline coefficients posterior distributions
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
sum_0s_single_50_gene <- ppc_stat(model_single_gene_50_score$y, yrep_model_gene_50, stat=function(x) sum(x==0)) + ggtitle('Sum of 0s') +
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
                                    fill = "Gene"
                                  )

# Clusterwise Models & plots ----------------------------------------------
#### Clusterwise fixed interaction (ns(Window)*Modules) ####
model_clusterwise_interaction <- readRDS("~/RNA-seq_training/scz_expression_bayes_regression_MCMC_rawcounts_clusterwise_interactions_TMM.rds")
yrep_clusterwise_interaction <- rstanarm::posterior_predict(model_clusterwise_interaction, draws = 500)
median_clusterwise_interaction <- ppc_stat(model_clusterwise_interaction$y, yrep_clusterwise_interaction, stat = "median", discrete = TRUE) + ggtitle('Median') + 
  theme(plot.title = element_text(hjust = 0.5, size = 12)) + theme(legend.position = "none")
sd_clusterwise_interaction <- ppc_stat(model_clusterwise_interaction$y, yrep_clusterwise_interaction, stat = "sd", discrete = TRUE) + ggtitle('Sd') + 
  theme(plot.title = element_text(hjust = 0.5, size = 12)) + theme(legend.position = "none")
sum_0s_clusterwise_interaction <- ppc_stat(model_clusterwise_interaction$y, yrep_clusterwise_interaction, stat=function(x) sum(x==0)) + ggtitle('Sum of 0s') +
  theme(plot.title = element_text(hjust = 0.5, size = 12))
