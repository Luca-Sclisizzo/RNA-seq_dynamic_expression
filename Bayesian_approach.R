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
})

# Parsing CLI arguments (normalization method and cores)
args <- commandArgs(trailingOnly = TRUE)
norm <- args[1]

normalization_methods <- c("TMM", "RLE", "upperquartile", "none")
if (! norm %in% normalization_methods) {
  stop(sprintf(
    "Unknown or missing normalization method: '%s'", norm
  ))
}
cores <- as.integer(args[2])
if (is.na(cores)) {
  cores <- parallel::detectCores()
}
#print(sprintf("Running in parallel over %d cores", cores))
#options(mc.cores = cores)
options(mc.cores = 1)

# Loading files
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
    expression = as.integer(round(.data$expression))
  ) %>%
  dplyr::filter(subject %in% sample_metadata$Braincode) # Keep only the EUR samples


##### Bayesian inference ##### 
# Bayesian inference
inference_result <- rstanarm::stan_glmer(
  expression ~ ns(Days, df = 4) + Sequencing.Site + Sex +
    (1 | area) +
    (1 | subject),
  data = int_RPKM_reshaped,  # logRPKM_reshaped,
  family = neg_binomial_2,
  #chains = min(, cores),
  #cores = 4,
  #adapt_delta = 0.8, 
  #control = list(max_treedepth =10),
  algorithm = "meanfield"
)

print('Done fitting the model, now saving the results...')
saveRDS(
  object = inference_result,
  file = paste0(
    "scz_expression_bayes_regression_meanfield_",
    norm,
    ".rds"
  )
)
