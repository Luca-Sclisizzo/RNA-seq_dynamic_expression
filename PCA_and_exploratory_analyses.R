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
  library("glmpca")
  library("patchwork")
  library("DESeq2")
})

normalization_methods <- c("TMM", "RLE", "upperquartile", "none")

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

##### Canonical PCA for QC #####
plots_PC12 <- list() # PC1 and PC2
plots_PC34 <- list() # PC3 and PC4

varimax <- FALSE
for (norm in normalization_methods){
  # group <- sub("\\..*$", "", colnames(counts)) # Braincode extraction
  dge <- edgeR::DGEList(counts = counts) #, group = group)
  keep <- edgeR::filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  
  # Normalize
  cat('\nNormalizing with',norm,'...')
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
    log = TRUE,
    prior.count = 1
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
  
  cat('\nStarting PCA analysis...')
  # PCA for different age samples
  df_wide <- int_RPKM_reshaped %>%
    mutate(sample = paste(subject, area, sep = ":")) %>%
    dplyr::select(sample, gene, expression) %>%
    pivot_wider(names_from = gene, values_from = expression)
  
  # Matrix for PCA
  mat <- df_wide %>%
    dplyr::select(-sample) %>%
    as.matrix()
  rownames(mat) <- df_wide$sample
  
  pca <- prcomp(mat,center = TRUE,scale. = TRUE)
  if(varimax == TRUE){ # Varimax rotation
    rot <- varimax(pca$rotation) # Loadings
    pca$rotation <- rot$loadings # Scores
    pca$x <- as.matrix(pca$x) %*% rot$rotmat
  }
  pca_df <- pca$x[,1:4] %>% # Selecting only the first 4 PCs
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    separate(sample, into = c("subject", "area"), sep = ":") %>%
    {
      if (varimax) {
        dplyr::rename(., PC1 = V1, PC2 = V2, PC3 = V3, PC4 = V4)
      } else {.}
    } %>%
    dplyr::select(subject, area, PC1, PC2, PC3, PC4) %>%
    left_join(sample_metadata %>% dplyr::select(Window, Braincode), by = c('subject' = 'Braincode'))
  
  cat('\nSaving the ggplot object...')
  plots_PC12[[norm]] <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Window)) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
    labs(x = "PC1", y = "PC2", title = norm) +
    theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  theme_minimal()
  
  plots_PC34[[norm]] <- ggplot(pca_df, aes(x = PC3, y = PC4, color = Window)) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
    labs(x = "PC3", y = "PC4", title = norm) +
    theme(plot.title = element_text(hjust = 0.5, size = 10)) +
    theme_minimal()
}


invisible(
  grid_plots_PC12 <- wrap_plots(plots_PC12, ncol = 2, nrow = 2) +
    plot_layout(guides = "collect") +
    plot_annotation(title = 'PC1 and PC2 - ROIs positions in the latent space colored by age') &
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5))
)
#print(grid_plots_PC12)


invisible(
  grid_plots_PC34 <- wrap_plots(plots_PC34, ncol = 2, nrow = 2) +
    plot_layout(guides = "collect") +
    plot_annotation(title = 'PC3 and PC4 - ROIs positions in the latent space colored by age') &
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5))
)
#print(grid_plots_PC34)

##### Generalized PCA for QC #####
generalized_PCA <- FALSE
if(generalized_PCA){
  plots <- list()
  for (norm in normalization_methods){
    # group <- sub("\\..*$", "", colnames(counts)) # Braincode extraction
    dge <- edgeR::DGEList(counts = counts) #, group = group)
    keep <- edgeR::filterByExpr(dge)
    dge <- dge[keep, , keep.lib.sizes=FALSE]
    
    # Normalize
    cat('\nNormalizing with',norm,'...')
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
    
    cat('\nStarting PCA analysis...')
    # PCA for different age samples
    df_wide <- int_RPKM_reshaped %>%
      mutate(sample = paste(subject, area, sep = ":")) %>%
      dplyr::select(sample, gene, expression) %>%
      pivot_wider(names_from = gene, values_from = expression)
    
    # Matrix for PCA
    mat <- df_wide %>%
      dplyr::select(-sample) %>%
      as.matrix()
    rownames(mat) <- df_wide$sample
    
    set.seed(40)
    gPCA <- glmpca(mat, L=2, fam = c('nb')) # generalized PCA
    
    pca_df <- gPCA$loadings %>%
      as.data.frame() %>%
      tibble::rownames_to_column("sample") %>%
      separate(sample, into = c("subject", "area"), sep = ":") %>%
      dplyr::select(subject, area, dim1, dim2) %>%
      left_join(sample_metadata %>% dplyr::select(Window, Braincode), by = c('subject' = 'Braincode'))
    
    cat('\nSaving the ggplot object...')
    plots[[norm]] <- ggplot(pca_df, aes(x = dim1, y = dim2, color = Window)) +
                          geom_point(size = 2) +
                          geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
                          geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
                          labs(x = "dim1", y = "dim2", title = norm) +
                          theme(plot.title = element_text(hjust = 0.5, size = 10))
                          theme_minimal()
  }
  invisible(
    grid_plots <- wrap_plots(plots, ncol = 2, nrow = 2) +
      plot_layout(guides = "collect") +
      plot_annotation(title = 'Generalized PCA - ROIs positions in the latent space colored by age') &
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5))
  )
  print(grid_plots)
}
