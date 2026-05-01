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

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)
gene_data <- getBM(
  attributes = c("ensembl_gene_id", "transcript_length"),
  mart = ensembl
)

scz_genes <- readxl::read_excel("SCZ_genes.xlsx")
# For each gene of interest, keep only the longest transcript
scz_genes <- scz_genes %>%
  dplyr::left_join(
    gene_data %>%
    by = c("GENE" = "ensembl_gene_id")
  ) %>%
  dplyr::group_by(GENE) %>%
  dplyr::slice_max(
    order_by = transcript_length,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

counts <- read.delim(
  "mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.count.txt",
  sep = "\t"
) %>%
  tidyr::separate_wider_delim(
    Geneid, delim = "|", names = c("ensembl_gene_id", "gene_name")
  ) %>%
  tibble::column_to_rownames("ensembl_gene_id") %>%
  select(-"gene_name")

group <- sub("\\..*$", "", colnames(counts))
dge <- edgeR::DGEList(counts = counts, group = group)
keep <- edgeR::filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

str(dge)
