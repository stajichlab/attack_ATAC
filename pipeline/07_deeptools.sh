#!/usr/bin/bash
#SBATCH -p short --mem 24gb -N 1 -n 24 --out logs/shared_sites.log

module load samtools
module load deeptools/3.4.3

CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
	CPU=1
fi


if [ $# -lt  7]; then
    echo "usage:sh 07_deeptools.sh <mping_sites_to_use>"
    echo "
    <output_folder_name>: this scripts will create an folder under current path;
    <peaks>: the output filtered peaks from choose_cutoff.sh, it should include three column, <chr> <start> <end> and it's 0 based;
    <density-file>: the ouput *.density file from the calculate_density step;
    <tag>: the output tag folder from prepareTag;
    <chrom_info>: a file include the chromosome name;
    <contig_prefix>: the name of the contig prefix. If there're no contig, set it as <no>;
    all the path should be the full path."
    exit -1; 
fi

#### mping_sites_to_use ####
tmpfoldername=$1


BAM=bam
MPING=data/mPing
TMP=tmp/$tmpfoldername

# genome lengths file
cut -f 1,2 $GENOME.fai > $GENOME.genome
COV=coverage/processed_shared
mkdir -p $COV $TMP


if [ ! -f $GENOME.fai ]; then
	module load samtools
	samtools faidx $GENOME
fi

for size in 250 500
do
	total=$(expr $size \* 2)
	for strain in HEG4 EG4 A123 A119
	do
		BEDFILE=$TMP/${strain}_shared.${total}_window.bed
		NEWBASE=${strain}_shared.${total}
        bamCoverage -b ${strain}.bam -o ${strain}.bw -bs=1 -p=max
        computeMatrix reference-point -S ${strain}.bw -R ${strain}_shared.${total}_window.bed --upstream 1000 --downstream 1000 -out ${strain}${total}_computeMatrix --sortUsing max --skipZeros -p=max --maxThreshold 100000
        plotHeatmap -m ${strain}${total}_computeMatrix -out ${strain}${total}_heatmap --colorList white,blue --sortUsing max
        plotProfile -m ${strain}${total}_computeMatrix -out ${strain}${total}_profile --perGroup --plotType=fill
	done
done
