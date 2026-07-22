# RNA-seq Training

This repository contains the materials and analyses developed for an RNA-seq training project.

The project focuses on transcriptomic data from the PsychENCODE consortium, which investigates gene expression patterns across different developmental stages of the human brain.

The dataset consists of post-mortem human brain biological samples.

## Data

Raw sequencing data are available through the PsychENCODE portal:

<http://development.psychencode.org/>

In this project, we focus on schizophrenia (SCZ)-associated genes identified by Trubetskoy et al. (GWAS):

<https://www.nature.com/articles/s41586-022-04434-5>

Gene mapping was performed using **hMAGMA**, integrating the following sources of information:

- positional mapping,
- expression quantitative trait loci (eQTL) information,
- Hi-C chromatin interaction data.

The resulting list of schizophrenia-associated genes is provided in:

`SCZ_genes.xlsx`

## Objectives

The main objectives of this training project are:

- Learn how to model RNA-seq data using Bayesian approaches.
- Perform quality control and preprocessing of transcriptomic data.
- Compare developmental expression trajectories across different normalization procedures and modelling strategies.
- Explore developmental patterns of gene expression in the human brain.

## Scripts Description

- **`PCA_and_exploratory_analyses.R`**\
  Performs exploratory data analysis and principal component analysis (PCA) to evaluate the impact of different normalization procedures on the structure of the dataset.

- **`Bayesian_approach.R`**\
  Implements hierarchical Bayesian models using the `rstanarm` package.

- **`plotting_results.R`**\
  Generates visualizations summarizing the results from hierarchical models fitted with `rstanarm`.

- **`Network_analysis.R`**\
  Performs weighted gene co-expression network analysis (WGCNA) using the `BioNERO` package.

- **`Bayesian_w_stan.R`**\
  Implements partial pooling hierarchical models using the `cmdstanr` interface.

- **`rstan_model.stan`**\
  Defines the Bayesian hierarchical model using Stan syntax. The model is compiled and executed from `Bayesian_w_stan.R`.

- **`PPC_rstan.R`**\
  Performs posterior predictive checks (PPCs) using the `posterior` package to summarize model diagnostics and generated quantities.

## Generated Data

The following files contain intermediate data generated during preprocessing and analysis:

- **`transcript_lengths.csv`**\
  Contains the lengths of all protein-coding transcripts corresponding to SCZ-associated genes. Transcript lengths were retrieved using the `biomaRt::useEnsembl()` function.

- **`gene_variability_score.csv`**\
  Contains gene-level quality scores, calculated as the ratio between inter-individual variance and intra-individual variance. Additional details are provided in the accompanying report.

- **`gene_network_membership.csv`**\
  Contains the membership of each gene to weighted co-expression network modules. These annotations are used as grouping structures in the partial pooling hierarchical models.

## Notes

This repository is intended for training and educational purposes.

## Running on HPC

### 1. Set up the R environment

Make sure that **R version 4.4.1** is available.

You can either:

- install R system-wide,
- load an appropriate HPC module,
- or create a dedicated conda environment:

``` bash
conda create -n my_r_env conda-forge::r=4.4.1
```

### 2. Initialize and restore the R virtual environment

Navigate to the repository directory, activate the conda environment, and start R:

``` bash
cd /path/to/here || exit 1
conda activate my_r_env
R
```

Inside R, restore the project dependencies using `renv`:

``` r
renv::restore()
```

If prompted with the following message:

``` text
It looks like you've called renv::restore() in a project that hasn't been activated yet.

How would you like to proceed?

1: Activate the project and use the project library.
2: Do not activate the project and use the current library paths.
3: Cancel and resolve the situation another way.

Selection:
```

Select option **1** to activate the project and use the project-specific library:

``` text
Selection: 1
```
