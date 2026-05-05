if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
BiocManager::install("edgeR")
Sys.unsetenv("BIOMART_HOST")
Sys.unsetenv("ENSEMBL_MART_HOST")

library("biomaRt")
library("edgeR")
library("dplyr")
library("tidyr")
library("readxl")
library("tibble")

# Loading files
transcript_lengths_file <- "transcript_lengths.csv"
scz_genes <- readxl::read_excel("SCZ_genes.xlsx")
sample_metadata <- read.csv("mRNA-seq_Sample metadata.csv", sep = ';', na.strings = c("NA", "")) %>%
  dplyr::filter(!is.na(Braincode) & Ethnicity == "European")


if (file.exists(transcript_lengths_file)) {
  print("Loading transcript lengths from file...")
  gene_data <- read.csv(transcript_lengths_file)
  
  gene_data <- gene_data %>%
    dplyr::rename(ensembl_gene_id = ensemble_gene_id) %>%
    dplyr::rename(transcript_length = trnascript_length) # Typo fix
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
  dplyr::left_join( # adding gene lenght
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
counts <- read.delim(
  "mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.count.txt",
  sep = "\t"
) %>%
  tidyr::separate_wider_delim(
    Geneid, delim = "|", names = c("ensembl_gene_id", "gene_name")
  ) %>%
  tibble::column_to_rownames("ensembl_gene_id") %>%
  select(-"gene_name")


##### Normalization factor for RNA-seq data ##### 
group <- sub("\\..*$", "", colnames(counts)) # Braincode extraction
dge <- edgeR::DGEList(counts = counts, group = group)
keep <- edgeR::filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Norm methods
#normalization_methods <- c("TMM","RLE","upperquartile","none") # possiamo pensare di parralelizzare questi jobs
normalization_methods <- c("TMM") # giusto per fare il primo run e vedere se funziona


for (norm in normalization_methods){ # let's see if the normalization method influences the transcriptomic trajectory
  dge <- edgeR::calcNormFactors(dge, method = norm)
  dge <- dge[rownames(dge) %in% scz_genes$GENE, ] # I have a drops of ~20 genes
  dge <- dge[, dge$samples$group %in% sample_metadata$Braincode] # Filtering only the EUR samples
  
  
  scz_genes <- scz_genes[ # reordering to avoid problems
    match(rownames(dge), scz_genes$GENE),
  ]
  
  logRPKM <- edgeR::rpkm(dge,
                        gene.length = scz_genes$transcript_length,
                        normalized.lib.size = TRUE, log = TRUE, prior.count = 1) # log and prior.count to avoid 0's problems

  
  # ANNOTATIONS: la funzione edgeR::rpkm() tiene conto dell'appartenza ai gruppi, lasciamo cosí?
  # O lasciamo che la gerarchia sia presa totalmente dal modello gerarchico downstream?
  # Mi chiedevo se mettendo group nella normalizzazione edgeR prendesse parte della "varianza gerarchica"
  # should we take into account also the zero inflation for the model?
  
  
}
