#!/usr/bin/env bash

cores="8"

for norm in "TMM" "RLE" "upperquartile" "none"; do
  sbatch  --time=24:00:00 --cpus-per-task="${cores}"  --mem=64G \
     Rscript ./preprocess_counts.R "${norm}" "${cores}"
done
