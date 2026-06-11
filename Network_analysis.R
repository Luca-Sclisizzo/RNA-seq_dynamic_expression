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
  library("BioNERO")
  library("viridis")
})
set.seed(123)
# The pipeline tutorial is available at https://bioconductor.org/packages//release/bioc/vignettes/BioNERO/inst/doc/vignette_01_GCN_inference.html#installation

# Data Loading & Preprocessing  -------------------------------------------
scz_genes <- readxl::read_excel("SCZ_genes.xlsx")
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
counts <- counts[rownames(counts) %in% scz_genes$GENE,] # filtering SCZ_genes

counts_preprocessed <- BioNERO::exp_preprocess( # preprocessing counts, this function is a wrapper of multiple preprocessing steps
  counts, min_exp = 10, variance_filter = FALSE, n = 2000
)

# Power Setting & Network Building -----------------------------------
sft <- SFT_fit(counts_preprocessed, net_type = "signed hybrid", cor_method = "pearson")
power <- sft$power # power to use the next step

net <- exp2gcn( # Network estimation
  counts_preprocessed, net_type = "signed hybrid", SFTpower = power, 
  cor_method = "pearson"
  )

#module_stability(counts_preprocessed, net, nRuns = 20)  + scale_fill_viridis_d(option = "cividis") # Assessing module stability
hubs <- get_hubs_gcn(counts_preprocessed, net) # Gene Hubs

# Plotting ----------------------------------------------------------------
net$genes_and_modules$Modules <- as.numeric(
  factor(net$genes_and_modules$Modules)
)
plot_genes_module <- plot_ngenes_per_module(net)




