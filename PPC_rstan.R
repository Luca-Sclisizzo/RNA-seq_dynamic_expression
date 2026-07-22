suppressPackageStartupMessages({
library('ggplot2')
library('dplyr')
library('splines')
library('bayesplot')
library('loo')
library('tidyr')
library('rstanarm')
library('patchwork')
library('cmdstanr')
library('posterior')
})

norm <- 'TMM'
par <- c(
  "alpha",
  as.vector(outer(1:9, 1:4, \(i, j) paste0("beta_module[", i, ",", j, "]")))
)

# RNA-seq preprocessing & data manipulation -------------------------------
transcript_lengths_file <- "transcript_lengths.csv"
gene_network_membership <- read.csv("gene_network_membership.csv")
scz_genes <- readxl::read_excel("SCZ_genes.xlsx")
sample_metadata <- read.csv(
  "mRNA-seq_Sample metadata.csv",
  sep = ';',
  na.strings = c("NA", "")
) %>%
  dplyr::filter(!is.na(Braincode) & Ethnicity == "European")

if (file.exists(transcript_lengths_file)) {
  print("Loading transcript lengths from file...")
  gene_data <- read.csv(transcript_lengths_file)
  
  gene_data <- gene_data %>%
    dplyr::rename(ensembl_gene_id = ensemble_gene_id)
} else {
  ensembl <- useEnsembl(
    biomart = "genes",
    dataset = "hsapiens_gene_ensembl"
  )
  gene_data <- getBM(
    attributes = c("ensembl_gene_id", "transcript_length"),
    mart = ensembl
  )
  print("Storing transcript lengths to file...")
  write.csv(
    gene_data %>%
      dplyr::filter(.data$ensembl_gene_id %in% scz_genes$GENE),
    transcript_lengths_file
  )
}

# For each gene of interest, keep only the longest transcript
scz_genes <- scz_genes %>%
  dplyr::left_join( # adding gene length
    gene_data,
    by = c("GENE" = "ensembl_gene_id")
  ) %>%
  dplyr::group_by(GENE) %>%
  dplyr::slice_max( # keep only the longest transcript
    order_by = transcript_length,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

# reading RNA-seq reads and preprocessing
# We only keep cortical areas (the ones ending in C, except CDC)
counts <- read.delim(
  "mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.count.txt",
  sep = "\t"
) %>%
  tidyr::separate_wider_delim(
    Geneid, delim = "|", names = c("ensembl_gene_id", "gene.name")
  ) %>%
  tibble::column_to_rownames("ensembl_gene_id") %>%
  dplyr::select(
    names(.)[
      sapply(strsplit(names(.), "\\."), \(x)
             endsWith(x[2], "C") && x[2] != "CBC"
      )
    ]
  )

# Normalization factor for RNA-seq data -----------------------------------
group <- sub("\\..*$", "", colnames(counts)) # Braincode extraction
dge <- edgeR::DGEList(counts = counts) #, group = group)
keep <- edgeR::filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Normalize
dge <- edgeR::calcNormFactors(dge, method = norm)
dge <- dge[rownames(dge) %in% scz_genes$GENE, ] # I have a drops of ~20 genes
# dge <- dge[, dge$samples$group %in% sample_metadata$Braincode] # Filtering only the EUR samples

scz_genes <- scz_genes[ # reordering to avoid problems
  match(rownames(dge), scz_genes$GENE),
]

# Raw Counts Manipulation -------------------------------------------------
int_RPKM <- dge$counts

normalization_factors <- dge$samples %>%
  tibble::rownames_to_column("sample") %>%
  mutate(offset = log(lib.size * norm.factors)) # this is the same than getOffset() function from edgeR

int_RPKM_reshaped <- as.data.frame(int_RPKM) %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "expression"
  ) %>%
  left_join(normalization_factors, by= 'sample') %>%
  tidyr::separate(sample, into = c("subject", "area"), sep = "\\.") %>%
  left_join(
    sample_metadata %>%
      dplyr::select("Braincode", "Days", "Sex", "Sequencing.Site", "Window"),
    by = c("subject" = "Braincode"),
  ) %>%
  left_join(
    scz_genes %>% dplyr::select(GENE, transcript_length),
    by = c('gene' = 'GENE')
  ) %>%
  mutate(
    transcript_length = transcript_length / 1000, # transcript length in kb
    expression = as.integer(round(.data$expression)),
    across(c(Sex, area, subject, gene,Sequencing.Site), as.factor)
  ) %>%
  dplyr::filter(subject %in% sample_metadata$Braincode) %>% # Keep only the EUR samples
    dplyr::left_join(gene_network_membership, by = c('gene' = 'Genes')) %>%
    subset(!is.na(Modules)) %>% # removing the genes that do not cluster to any group
    mutate(Modules = as.factor(Modules)) %>% # set it as a factor
  mutate(
    gene_idx    = as.integer(factor(gene)),
    module_idx  = as.integer(factor(Modules))
  )

# gene Mapping
gene_module_map <- int_RPKM_reshaped %>%
  distinct(gene_idx, module_idx) %>%
  arrange(gene_idx) %>%
  pull(module_idx)


# rstan data list creation -------------------------------------------

B <- ns(int_RPKM_reshaped$Window, df = 4) # building the spline
B <- scale(B)
subject <- as.integer(factor(int_RPKM_reshaped$subject)) # Creating J groups 1..J-th

stan_data <- list(
  N = nrow(int_RPKM_reshaped),
  K = ncol(B),
  S = length(unique(subject)),
  M = length(unique(int_RPKM_reshaped$Modules)),
  G = length(unique(int_RPKM_reshaped$gene)),
  
  module      = int_RPKM_reshaped$module_idx,
  gene        = int_RPKM_reshaped$gene_idx,
  gene_module = gene_module_map,
  age         = int_RPKM_reshaped$Window,
  subject     = subject,
  log_offset = int_RPKM_reshaped$offset,
  
  B = B,
  y = int_RPKM_reshaped$expression
)

# PPC graphic & Saving -------------------------------
results_list <- list()
model_clusterwise_interaction_w_variability <- readRDS('./scz_expression_bayes_regression_rstan_clusterwise_gene_variability_TMM.rds')
print('Loading the model for y_rep generation...')
model <- cmdstan_model("rstan_model.stan") # model with y_rep generation

convergence_checks <- model_clusterwise_interaction_w_variability$summary(variables = par) %>%
  dplyr::filter(rhat > 1.01 | ess_bulk < 400)

# Neff ratio calculation and plotting
draws_arr <- model_clusterwise_interaction_w_variability$draws(variables = par)
n_draws <- niterations(draws_arr) * nchains(draws_arr)
ess_summary <- summarise_draws(draws_arr, ess_bulk)
neff_ratios <- ess_summary$ess_bulk / n_draws
names(neff_ratios) <- ess_summary$variable

neff_plot <- mcmc_neff(neff_ratios) +
  ggtitle("Neff Ratio Plot") + theme(plot.title = element_text(hjust = 0.5, size = 12))

# MCMC dens overlay plot
dens_plot <- mcmc_dens_overlay(
  model_clusterwise_interaction_w_variability$draws(format = "df"),
  pars = par[1:6]
) +
  ggtitle("Density per Chain") + theme(plot.title = element_text(hjust = 0.5, size = 12))

# y_rep generation and PPC plots
draws_df <- model_clusterwise_interaction_w_variability$draws(
   format = "draws_matrix",
 )[1:100, ]  # first 100 draws

gq_fit <- model$generate_quantities(
   fitted_params = draws_df,
   data = stan_data
)
y_rep <- gq_fit$draws("y_rep", format = "matrix")  # [S x N]

mean_ppc <- ppc_stat(stan_data$y, y_rep, stat = "median", discrete = TRUE) + ggtitle('Median') +
  theme(plot.title = element_text(hjust = 0.5, size = 12)) + theme(legend.position = "none")
sd_ppc <- ppc_stat(stan_data$y, y_rep, stat = "sd", discrete = TRUE) + ggtitle('Sd') +
  theme(plot.title = element_text(hjust = 0.5, size = 12)) + theme(legend.position = "none")
sum_0s <- ppc_stat(stan_data$y, y_rep, stat=function(x) sum(x==0)) + ggtitle('Sum of 0s') +
  theme(plot.title = element_text(hjust = 0.5, size = 12))
PPC_list <- list(mean_ppc, sd_ppc, sum_0s)

PPC_patchwork <- patchwork::wrap_plots(PPC_list) +
     plot_annotation(title = 'PPC in partial pooling with gene-wise intercepts',
     theme = theme(plot.title = element_text(hjust = 0.5), 
     plot.caption = element_text(hjust = 0)
))

# Saving results
results_list$PPC <- PPC_patchwork
results_list$neff_plot <- neff_plot
results_list$dens_plot <- dens_plot
results_list$convergence_checks <- convergence_checks
saveRDS(results_list, './results_list.rds')
