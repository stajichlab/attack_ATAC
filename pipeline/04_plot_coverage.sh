#!/usr/bin/bash
#SBATCH -p short --out logs/plot_coverage.log
mkdir -p plots
Rscript plot_coverage.R
