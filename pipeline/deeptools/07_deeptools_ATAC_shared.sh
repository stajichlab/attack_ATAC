#!/usr/bin/bash -l
#SBATCH -p short --mem 24gb -N 1 -n 24 --out logs/deeptools_ATAC_shared.log

module load samtools
module unload miniconda2
module unload miniconda3
module load anaconda3
module load deeptools/3.4.3
source activate deeptools

CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
	CPU=1
fi


if [ $# -lt  3 ]; then
    echo "usage:sh 07_deeptools.sh <mping_sites_to_use> <peaks> <density-file> <tags> <chrom_info> <contig_prefix>"
    echo "
    <datatype>: ATAC
    <tmpfoldername>: 03_processed_empty or 04_background_sites or 05_processed_shared
    <processed>: processed_background or processed_empty or processed_shared

    all the path should be the full path."
    exit -1; 
fi

#### mping_sites_to_use ####
datatype=$1
tmpfoldername=$2
processed=$3
TMP=tmp/$tmpfoldername
BAM=bam
MPING=data/mPing

COV=coverage/$processed/deeptools
# make folder if it doesn't exist 
mkdir -p $COV 

    
for size in 250 500
    do
        total=$(expr $size \* 2)
        for strain in HEG4 EG4 A123 A119
        do
            BEDFILE=$TMP
            NEWBASE=${strain}_shared.
            # I think this currently runs but it's ridiculously slow for some reason
            bamCoverage -b $BAM/$datatype/${strain}.bam -o $COV/$datatype${strain}.bw -bs=1 -p=max --normalizeUsing RPKM

            # Shared
            computeMatrix reference-point -S $COV/$datatype${strain}.bw -R tmp/05_processed_shared/${strain}_shared_sites.bed   --upstream 2000 --downstream 2000 -out $COV/$datatypeCMref${strain}_shared --binSize 50 --sortUsing max --skipZeros -p=max --maxThreshold 100000
            plotHeatmap -m $COV/$datatypeCMref${strain}_shared -out $COV/${datatype}CMref${strain}_sharedheatmap --colorList white,blue --sortUsing max --plotTitle ${datatype}CMref${strain}_shared --refPointLabel mPing
            plotHeatmap -m $COV/$datatypeCMref${strain}_shared -out $COV/${datatype}CMref${strain}_sharedheatmapkmeans --outFileSortedRegions $COV/${datatype}CMref${strain}_sharedheatmapkmeans.stats --colorList white,blue --sortUsing max --plotTitle ${datatype}CMref${strain}_shared --kmeans 6 --refPointLabel mPing    
            plotProfile -m $COV/$datatypeCMref${strain}_shared -out $COV/${datatype}CMref${strain}_sharedprofile --perGroup --plotType=se --plotTitle ${datatype}CMref${strain}_shared  --refPointLabel mPing
            plotProfile -m $COV/$datatypeCMref${strain}_shared -out $COV/${datatype}CMref${strain}_sharedprofilekmeans --outFileSortedRegions $COV/${datatype}CMref${strain}_sharedprofilekmeans.stats --perGroup --plotType=se --plotTitle ${datatype}CMref${strain}_shared --kmeans 6 --refPointLabel mPing

   done
done