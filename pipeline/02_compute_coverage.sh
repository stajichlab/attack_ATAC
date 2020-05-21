#!/usr/bin/bash
#SBATCH -p short --mem 64gb -N 1 -n 24 --out logs/coverage.log

module load bedtools
module unload perl
module load parallel
module load mosdepth

CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
	CPU=1
fi
CPUS=$CPU

BAM=bam
GENOME=genome/Oryza_sativa.IRGSP-1.0.30.dna.toplevel.fa
COV=coverage/processed_windows
mkdir -p $COV
# could do this with parallel with 2 arguments but that's prob even more confusing...
for WINDOW in 10000 5000 1000 500
do
    parallel --rpl '{///} $Global::use{"File::Basename"} ||= eval "use File::Basename; 1;"; $_ = basename(dirname($_));' \
		 -j $CPUS mosdepth -f $GENOME -x -n --by $WINDOW $COV/{///}_{/.}_w${WINDOW} {} ::: $(ls $BAM/*/*.bam)
done

mkdir coverage/processed_windows/tmp
mv coverage/processed_windows/*bed.gz.csi coverage/processed_windows/tmp
mv coverage/processed_windows/*global.dist.txt coverage/processed_windows/tmp
mv coverage/processed_windows/*region.dist.txt coverage/processed_windows/tmp

mkdir coverage/processed_windows/ATAC
mkdir coverage/processed_windows/H3K27me3
mkdir coverage/processed_windows/H3K36me3
mkdir coverage/processed_windows/H3K56ac

mv coverage/processed_windows/ATAC* coverage/processed_windows/ATAC
mv coverage/processed_windows/H3K27me3* coverage/processed_windows/H3K27me3
mv coverage/processed_windows/H3K36me3* coverage/processed_windows/H3K36me3
mv coverage/processed_windows/H3K56ac* coverage/processed_windows/H3K56ac