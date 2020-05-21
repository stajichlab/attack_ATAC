#!/usr/bin/bash
#SBATCH -p short --out logs/plot_coverage.log

Rscript plot_coverage.R
