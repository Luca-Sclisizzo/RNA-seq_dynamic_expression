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

norm <- 'TMM'
normalization_methods <- c("TMM", "RLE", "upperquartile", "none")
if (! norm %in% normalization_methods) {
  stop(sprintf(
    "Unknown or missing normalization method: '%s'", norm
  ))
}
args <- commandArgs(trailingOnly = TRUE)
cores <- as.integer(args[1])
if (is.na(cores)) {
  cores <- parallel::detectCores()
}
print(sprintf("Running in parallel over %d cores", cores))
options(mc.cores = cores)

# RNA-seq preprocessing & data manipulation -------------------------------
transcript_lengths_file <- "transcript_lengths.csv"
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
# group <- sub("\\..*$", "", colnames(counts)) # Braincode extraction
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
    Sequencing.Site = as.factor(Sequencing.Site)
  ) %>%
  dplyr::filter(subject %in% sample_metadata$Braincode) # Keep only the EUR samples

# Bayesian Inference with rstan -------------------------------------------
print("Fitting the model with rstan...\n")
print("Compiling the model...")
model <- cmdstan_model("rstan_model.stan")

B <- ns(int_RPKM_reshaped$Window, df = 4) # building the spline
B <- scale(B)
subject <- as.integer(factor(int_RPKM_reshaped$subject)) # Creating J groups 1..J-th

stan_data <- list( # shaping the df as a list of parameters
  N = nrow(int_RPKM_reshaped),
  K = ncol(B),
  S = length(unique(subject)),
  
  age = int_RPKM_reshaped$Window,
  subject = subject, # Creating J groups 1..J-th
  
  B = B, # Spline
  y = int_RPKM_reshaped$expression # RNA-seq
)

fit <- model$sample( # Fitting the model 
  data = stan_data,
  iter_sampling = 100,
  iter_warmup = 100,
  chains = 1,
  parallel_chains = 4,
  init = 0.5
)

# Saving ------------------------------------------------------------------
print('Done fitting the model, now saving the results...')
fit$save_object(
  paste0(
    "scz_expression_bayes_regression_rstan_",
    norm,
    ".rds"
  ))
