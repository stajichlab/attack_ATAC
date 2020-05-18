#!/usr/bin/bash
#SBATCH -p short -N 1 -n 2 --mem 24gb

bash pipeline/00_download.sh
bash pipeline/01_make_windows.sh
sbatch pipeline/02_compute_coverage.sh
sbatch pipeline/03_empty_sites.sh

