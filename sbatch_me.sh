#!/usr/bin/env bash

for norm in "TMM" "RLE" "upperquartile" "none"; do
  sbatch  --time=24:00:00 --cpus-per-task=8  --mem=64G \
     Rscript ./preprocess_counts.R "$norm"
done
