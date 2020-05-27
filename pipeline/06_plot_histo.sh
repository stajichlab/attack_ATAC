#!/usr/bin/bash
#SBATCH -p short --mem 4gb -N 1 -n 8 --out logs/plot_R.log

Rscript Rscripts/plot_coverage.R
