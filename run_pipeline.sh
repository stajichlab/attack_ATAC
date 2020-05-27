#!/usr/bin/bash
#SBATCH -p short -N 1 -n 2 --mem 16gb 

bash pipeline/00_download.sh
bash pipeline/01_make_windows.sh
sbatch pipeline/02_compute_coverage.sh
sbatch pipeline/03_empty_sites.sh
sbatch pipeline/04_background_sites.sh  
 
# when all jobs are done then you can run 06_plot_coverage.sh
# we can make this happen using slurm dependencies but I'll leave that out for now
