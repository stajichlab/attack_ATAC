#!/usr/bin/bash
#SBATCH -p short -N 1 -n 2 --mem 16gb --out logs/deeptools.log

sbatch ./pipeline/deeptools/07_deeptools_ATAC_shared.sh ATAC 05_processed_shared processed_shared
sbatch ./pipeline/deeptools/07_deeptools_ChIP_shared.sh H3K27me3 05_processed_shared processed_shared
sbatch ./pipeline/deeptools/07_deeptools_ChIP_shared.sh H3K36me3 05_processed_shared processed_shared
sbatch ./pipeline/deeptools/07_deeptools_ChIP_shared.sh H3K56ac 05_processed_shared processed_shared

sbatch ./pipeline/deeptools/07_deeptools_ATAC_empty.sh ATAC 03_processed_empty processed_empty
sbatch ./pipeline/deeptools/07_deeptools_ChIP_empty.sh H3K27me3 03_processed_empty processed_empty
sbatch ./pipeline/deeptools/07_deeptools_ChIP_empty.sh H3K36me3 03_processed_empty processed_empty
sbatch ./pipeline/deeptools/07_deeptools_ChIP_empty.sh H3K56ac 03_processed_empty processed_empty

sbatch ./pipeline/deeptools/07_deeptools_ATAC_background.sh ATAC 04_background_sites processed_background
sbatch ./pipeline/deeptools/07_deeptools_ChIP_background.sh H3K27me3 04_background_sites processed_background
sbatch ./pipeline/deeptools/07_deeptools_ChIP_background.sh H3K36me3 04_background_sites processed_background
sbatch ./pipeline/deeptools/07_deeptools_ChIP_background.sh H3K56ac 04_background_sites processed_background
