#!/usr/bin/bash -l
#SBATCH -p short --mem 24gb -N 1 -n 24 --out logs/deeptoolsATAC_shared.log

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
TMP = tmp/$tmpfoldername
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
            NEWBASE=${strain}_shared.${total}
            # I think this currently runs but it's ridiculously slow for some reason
            #bamCoverage -b $BAM/$datatype/${strain}.bam -o $COV/$datatype${strain}.bw -bs=1 -p=max

            # Background
            #computeMatrix reference-point -S $COV/$datatype${strain}.bw -R tmp/04_background_sites/Simulate0001_${strain}_background_sites.bed   --upstream 2000 --downstream 2000 -out $COV/$datatypeCMref${strain}_background --binSize 50 --sortUsing max --skipZeros -p=max --maxThreshold 100000
            #plotHeatmap -m $COV/$datatypeCMref${strain}_background -out $COV/${datatype}CMref${strain}_backgroundheatmap --colorList white,blue --sortUsing max
            #plotProfile -m $COV/$datatypeCMref${strain}_background -out $COV/${datatype}CMref${strain}_backgroundprofile --perGroup --plotType=fill

            # Shared
            computeMatrix reference-point -S $COV/$datatype${strain}.bw -R tmp/05_processed_shared/${strain}_shared_sites.bed   --upstream 2000 --downstream 2000 -out $COV/$datatypeCMref${strain}_shared --binSize 50 --sortUsing max --skipZeros -p=max --maxThreshold 100000
            plotHeatmap -m $COV/$datatypeCMref${strain}_shared -out $COV/${datatype}CMref${strain}_sharedheatmap --colorList white,blue --sortUsing max
            plotProfile -m $COV/$datatypeCMref${strain}_shared -out $COV/${datatype}CMref${strain}_sharedprofile --perGroup --plotType=fill
            
   done
done