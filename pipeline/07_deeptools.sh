#!/usr/bin/bash
#SBATCH -p short --mem 24gb -N 1 -n 24 --out logs/deeptools.log

module load samtools
module load deeptools/3.4.3
module unload python/2.7.5
module load miniconda3
source activate deeptools

CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
	CPU=1
fi


if [ $# -lt  3]; then
    echo "usage:sh 07_deeptools.sh <mping_sites_to_use> <peaks> <density-file> <tags> <chrom_info> <contig_prefix>"
    echo "
    <Data_analyzed>: ATAC or H3K27me3 or H3K36me3 or H3K56ac
    <output_folder_name>: 03_processed_empty or 04_background_sites or 05_processed_shared
    <processedfolder>: processed_background or processed_empty or processed_shared

    all the path should be the full path."
    exit -1; 
fi

#### mping_sites_to_use ####
datatype=$1
tmpfoldername=$2
processed=$3

BAM=bam
MPING=data/mPing
TMP=tmp/$tmpfoldername
COV=coverage/$processed/deeptools
# make folder if it doesn't exist 
mkdir -p $COV $TMP


for size in 250 500
do
	total=$(expr $size \* 2)
	for strain in HEG4 EG4 A123 A119
	do
		BEDFILE=$TMP/$tmpfoldername
		NEWBASE=${strain}_shared.${total}
        # I think this currently runs but it's ridiculously slow for some reason
        bamCoverage -b $BAM/$datatype/${strain}.bam -o $COV/${strain}.bw -bs=1 -p=max
        computeMatrix reference-point -S $COV/${strain}.bw -R $BEDFILE --upstream 2000 --downstream 2000 -out $COV/CMref${strain}${total} --binSize 50 --sortUsing max --skipZeros -p=max --maxThreshold 100000
        plotHeatmap -m $COV/CMref${strain}${total} -out $COV/CMref${strain}${total}heatmap --colorList white,blue --sortUsing max
        plotProfile -m $COV/CMref${strain}${total} -out $COV/CMref${strain}${total}profile --perGroup --plotType=fill
	done
done