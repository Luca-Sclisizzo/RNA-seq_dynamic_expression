# RNA-seq Training

This repository contains materials and analyses for an RNA-seq training project.

We analyse transcriptomic data from the PsychENCODE consortium, which investigates gene expression across different developmental stages of the human brain.

Data comes from post-mortem biological samples.

## Data
Raw sequencing reads are available at:
http://development.psychencode.org/#

We will analyze SCZ genes from Trubetskoy et al. GWAS (https://www.nature.com/articles/s41586-022-04434-5).

hMAGMA analysis is performed to map the SNPs into genes, using both positional mapping, eQTLs and HiC interactions.
The list of Schizophrenia genes is in the SCZ_genes.xlsx file.

## Objectives
- Learn to model RNA-seq data with Bayesian approach
- Perform quality control and preprocessing
- Compare trajectories across normalization methods and modelling solutions 
- Explore developmental transcriptomic patterns

## Scripts description
- PCA_and_exploratory_analyses.R wants to perform exploratory analyses and PCA to better understand the outcome of different normalization processes
- Frequentist_approach.R performs hierarchical models in a traditional way, checking convergence problems and so on. This is necessary due to computational burden of bayesian samplings
- Bayesian_approach.R performs bayesian hierarchical models
- plotting_results.R plots the results of hierarchical models
## Notes
This repository is intended for training and educational purposes.

## How to run on HPC
* Make sure you have `Rv4.4.1`
    * Install it system-wide
    * Load a corresponding HPC module
    * Install it inside a conda environment with `conda create -n my_r_env
      conda-forge::r=4.4.1`
* Initialize and restore the R virtual environment
```bash
cd /path/to/here || exit 1
conda activate my_r_env
R
> renv::restore()
> It looks like you've called renv::restore() in a project that hasn't been activated yet.
How would you like to proceed?

1: Activate the project and use the project library.
2: Do not activate the project and use the current library paths.
3: Cancel and resolve the situation another way.

Selection: 1
```
* Run the [`sbatch_me.sh`](./sbatch_me.sh) script, after adjusting it to reflect
  the path you cloned this repo in and the absolute path of your Rscript
  interpreter, e.g.:
  ```bash
    /path/to/miniconda3/envs/R4.4/bin/Rscript \
        /path/to/RNA-seq_training/preprocess_counts.R \
        "$norm"
  ```
