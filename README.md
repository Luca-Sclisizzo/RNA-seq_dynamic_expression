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

## Notes
This repository is intended for training and educational purposes.