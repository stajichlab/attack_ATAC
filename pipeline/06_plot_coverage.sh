#!/usr/bin/bash
#SBATCH -p short --out logs/plot_coverage.log
#mkdir -p plots
Rscript Rscripts/plot_coverage.R
Rscript Rscripts/KS_compare_sites.R
