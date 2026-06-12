if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- c("biomaRt", "edgeR")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p)
  }
}
Sys.unsetenv("BIOMART_HOST")
Sys.unsetenv("ENSEMBL_MART_HOST")

# Setting a single thread for all the computations to avoid issues with parallelization in rstanarm
#Sys.setenv(OMP_NUM_THREADS = 1)
#Sys.setenv(OPENBLAS_NUM_THREADS = 1)
#Sys.setenv(MKL_NUM_THREADS = 1)


suppressPackageStartupMessages({
  library("biomaRt")
  library("dplyr")
  library("edgeR")
  library("readxl")
  library("rstanarm")
  library("splines")
  library("tibble")
  library("tidyr")
  library("lme4")
  library("ggplot2")
  library("cmdstanr")
})

# Parsing CLI arguments (normalization method and cores)
args <- commandArgs(trailingOnly = TRUE)
norm <- args[1]
n_gene_to_fit <- as.integer(args[3]) # this will select the n-th gene based on the gene_variability score ranking
normalization_methods <- c("TMM", "RLE", "upperquartile", "none")
if (! norm %in% normalization_methods) {
  stop(sprintf(
    "Unknown or missing normalization method: '%s'", norm
  ))
}
clusterwise_model <- TRUE
raw_counts <- TRUE # this is a switch between raw counts and pseudocounts
if (raw_counts) {
  warning("raw_counts = TRUE: analyses will use raw counts, library size and library composition are used as offset in the NB")
  if(clusterwise_model){
    warning("clusterwise_model = TRUE: model will use a partial pooling across clusters")
  }
} else {
  warning("raw_counts = FALSE: analyses will use pseudocounts, no offsets added to the NB")
  if(clusterwise_model == TRUE){
    stop("It is not possible to run a clusterwise model with pseudocounts.\nSet raw_counts = TRUE")
  }
}
cores <- as.integer(args[2])
if (is.na(cores)) {
  cores <- parallel::detectCores()
}
print(sprintf("Running in parallel over %d cores", cores))
options(mc.cores = cores)
#options(mc.cores = 1)

# Loading files
transcript_lengths_file <- "transcript_lengths.csv"
gene_variability <- read.csv("gene_variability_score.csv") # this file contains the gene_wise score of var_between_window / var_within_window
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

##### Normalization factor for RNA-seq data ##### 
group <- sub("\\..*$", "", colnames(counts)) # Braincode extraction
dge <- edgeR::DGEList(counts = counts, group = group)
keep <- edgeR::filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Normalize
dge <- edgeR::calcNormFactors(dge, method = norm)
dge <- dge[rownames(dge) %in% scz_genes$GENE, ] # I have a drops of ~20 genes
# dge <- dge[, dge$samples$group %in% sample_metadata$Braincode] # Filtering only the EUR samples

scz_genes <- scz_genes[ # reordering to avoid problems
  match(rownames(dge), scz_genes$GENE),
]
##### Bayesian inference - pseudocounts ##### 
if(raw_counts == FALSE){
  # If we want to use integers, we set prior.count to 0 and set log = FALSE
  # Problem with this approach: these normalized pseudo-counts are not generated
  #   by a true negative binomial model as the raw counts are.
  int_RPKM <- edgeR::rpkm(
    dge,
    gene.length = scz_genes$transcript_length,
    normalized.lib.size = TRUE,
    log = FALSE,
    prior.count = 0
  )
  int_RPKM_reshaped <- as.data.frame(int_RPKM) %>%
    rownames_to_column("gene") %>%
    pivot_longer(
      cols = -gene,
      names_to = "sample",
      values_to = "expression"
    ) %>%
    tidyr::separate(sample, into = c("subject", "area"), sep = "\\.") %>%
    left_join(
      sample_metadata %>%
        dplyr::select("Braincode", "Days", "Sex", "Sequencing.Site", "Window"),
      by = c("subject" = "Braincode"),
    ) %>%
    dplyr::mutate(
      expression = as.integer(round(.data$expression)),
      across(c(Sex, area, subject, gene, Sequencing.Site), as.factor)
    ) %>%
    dplyr::filter(subject %in% sample_metadata$Braincode) # Keep only the EUR samples
  
  # data QC
  # zero_rate <- int_RPKM_reshaped %>%
  #   group_by(gene) %>%
  #   summarise(prop_zero = mean(expression == 0))
  # 
  # counts_for_lowly_expressed <- int_RPKM_reshaped %>%
  #   subset(gene %in% zero_rate[zero_rate$prop_zero >= 0.42,]$gene) # I chose .42 because the max expression was <8 and the 3rd quartile was 1
  # 
  # int_RPKM_reshaped <- int_RPKM_reshaped %>%
  #   subset(!gene %in% counts_for_lowly_expressed$gene) # keep only the trustworthy genes
  # hist(log(int_RPKM_reshaped$expression + 1))
  
  # Bayesian inference
  inference_result <- rstanarm::stan_glmer(
    expression ~ ns(Window, df = 4) + Sex + Sequencing.Site + area + # fixed effects
      (1 | subject), # random effects
    data = int_RPKM_reshaped,
    family = neg_binomial_2,
    chains = min(4, cores),
    cores = cores,
    adapt_delta = 0.8, # This is because rstanarm is more conservative than rstan, this is the rstan default value
    control = list(max_treedepth =10), # This is because rstanarm is more conservative than rstan, this is the rstan default value
    algorithm = "sampling",
    iter = 3000,
    warmup = 2000
  )
  print('Done fitting the model, now saving the results...')
  saveRDS(
    object = inference_result,
    file = paste0(
      "scz_expression_bayes_regression_MCMC_gene_by_window_",
      norm,
      ".rds"
    )
  )
}

##### Bayesian inference - Raw Counts #####
# This section wants to try to use the raw counts and normalize in the model instead of using pseudocounts
# The idea is to use the library size and the gene length as an exposure measure
if(raw_counts == TRUE){ # only if raw counts are selected
  int_RPKM <- dge$counts

  normalization_factors <- dge$samples %>%
    rownames_to_column("sample") %>%
    mutate(offset = log(lib.size * norm.factors)) # this is the same than getOffset() function from edgeR
  
  int_RPKM_reshaped <- as.data.frame(int_RPKM) %>%
    rownames_to_column("gene") %>%
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
      scz_genes %>% select(GENE, transcript_length),
      by = c('gene' = 'GENE')
    ) %>%
    mutate(
      transcript_length = transcript_length / 1000, # transcript length in kb
      expression = as.integer(round(.data$expression)),
      across(c(Sex, area, subject, gene,Sequencing.Site), as.factor)
      ) %>%
    dplyr::filter(subject %in% sample_metadata$Braincode) # Keep only the EUR samples

  if (exists("n_gene_to_fit") && !is.null(n_gene_to_fit) && !is.na(n_gene_to_fit) && clusterwise_model == FALSE) {
    warning(sprintf(
      "Selected the %d-th gene based on the gene_variability score ranking from the gene_variability score.",
      n_gene_to_fit
    ))
    warning("A single gene model will be fitted!")
    int_RPKM_reshaped <- int_RPKM_reshaped %>% 
      dplyr::filter(gene %in% gene_variability$gene[n_gene_to_fit])
  } else if (exists("n_gene_to_fit") && !is.null(n_gene_to_fit) && !is.na(n_gene_to_fit) && clusterwise_model == TRUE) {
    stop("clusterwise_model = TRUE & n_gene_to_fit specified. Impossible to perform a clusterwise model on a single gene.")
  }
  if(!exists("n_gene_to_fit") && clusterwise_model == TRUE){
    warning("Clusterwise model will be perfomed. Adding cluster informations...")
    int_RPKM_reshaped <- int_RPKM_reshaped %>% 
      dplyr::left_join(gene_network_membership, by = c('gene' = 'Genes')) %>%
      subset(!is.na(Modules)) %>% # removing the genes that do not cluster to any group
      mutate(Modules = as.factor(Modules)) # set it as a factor
  }
  
  # Model using raw counts and normalization as an offset
  # The offset will be the edgeR offset: log(lib.size * norm_factor), it can be retreived with edgeR::getOffset(dge) function
  # Bayesian inference
  if(clusterwise_model == FALSE){
    inference_result <- rstanarm::stan_glmer(
      expression ~ ns(Window, df = 4) + Sex + Sequencing.Site + area + # fixed effects
        (1 | subject), # random effects
      data = int_RPKM_reshaped,
      family = neg_binomial_2,
      offset = offset,
      chains = min(4, cores),
      cores = cores,
      adapt_delta = 0.8, # rstanarm is more conservative than rstan, this is the rstan default value
      control = list(max_treedepth =10), # rstanarm is more conservative than rstan, this is the rstan default value
      algorithm = "sampling"#,
      #iter = 1000,
      #warmup = 2000
    )
    print('Done fitting the model, now saving the results...')
    saveRDS(
      object = inference_result,
      file = paste0(
        "scz_expression_bayes_regression_MCMC_rawcounts_single_gene_",
        n_gene_to_fit,
        "_th_gene_",
        norm,
        ".rds"
      )
    )
  } else { # Clusterwise bayesian inference
    inference_result <- rstanarm::stan_glmer(
      expression ~ ns(Window, df = 4) + Sex + Sequencing.Site + area + # fixed effects
        (1 | subject) + (1 | Modules), # random effects
      data = int_RPKM_reshaped,
      family = neg_binomial_2,
      offset = offset,
      chains = min(4, cores),
      cores = cores,
      adapt_delta = 0.8, # rstanarm is more conservative than rstan, this is the rstan default value
      control = list(max_treedepth =10), # rstanarm is more conservative than rstan, this is the rstan default value
      algorithm = "sampling"#,
      #iter = 1000,
      #warmup = 2000
    )
    print('Done fitting the model, now saving the results...')
    saveRDS(
      object = inference_result,
      file = paste0(
        "scz_expression_bayes_regression_MCMC_rawcounts_clusterwise_",
        norm,
        ".rds"
      )
    )
  }
}
print('Done fitting and saving, job completed!')
