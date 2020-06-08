#!/usr/bin/bash -l
#SBATCH -p short --mem 24gb -N 1 -n 24 --out logs/deeptools.log

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
    <datatype>: ATAC or H3K27me3 or H3K36me3 or H3K56ac
    <tmpfoldername>: 03_processed_empty or 04_background_sites or 05_processed_shared
    <processed>: processed_background or processed_empty or processed_shared

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
mkdir -p $COV 
mkdir -p $TMP

if [ $datatype = "ATAC" ]; then
    for size in 250 500
    do
        total=$(expr $size \* 2)
        for strain in HEG4 EG4 A123 A119
        do
            BEDFILE=$TMP
            NEWBASE=${strain}_shared.${total}
            # I think this currently runs but it's ridiculously slow for some reason
            bamCoverage -b $BAM/$datatype/${strain}.bam -o $COV/$datatype${strain}.bw -bs=1 -p=max

            if [ $tmpfoldername = "04_background_sites" ]; then
           		# Background
           		computeMatrix reference-point -S $COV/$datatype${strain}.bw -R tmp/04_background_sites/Simulate0001_${strain}_background_sites.bed   --upstream 2000 --downstream 2000 -out $COV/$datatypeCMref${strain}${total}_background --binSize 50 --sortUsing max --skipZeros -p=max --maxThreshold 100000
            	plotHeatmap -m $COV/$datatypeCMref${strain}${total}_background -out $COV/${datatype}CMref${strain}${total}_backgroundheatmap --colorList white,blue --sortUsing max
            	plotProfile -m $COV/$datatypeCMref${strain}${total}_background -out $COV/${datatype}CMref${strain}${total}_backgroundprofile --perGroup --plotType=fill
			elif [ $tmpfoldername = "05_processed_shared"  ];
            	# Shared
            	computeMatrix reference-point -S $COV/$datatype${strain}.bw -R tmp/05_processed_shared/${strain}_shared_sites.bed   --upstream 2000 --downstream 2000 -out $COV/$datatypeCMref${strain}${total}_shared --binSize 50 --sortUsing max --skipZeros -p=max --maxThreshold 100000
            	plotHeatmap -m $COV/$datatypeCMref${strain}${total}_shared -out $COV/${datatype}CMref${strain}${total}_sharedheatmap --colorList white,blue --sortUsing max
            	plotProfile -m $COV/$datatypeCMref${strain}${total}_shared -out $COV/${datatype}CMref${strain}${total}_sharedprofile --perGroup --plotType=fill
			elif [ $tmpfoldername = "03_processed_empty" ];
            	# Empty
            	computeMatrix reference-point -S $COV/$datatype${strain}.bw -R tmp/03_processed_empty/${strain}_empty_sites.bed   --upstream 2000 --downstream 2000 -out $COV/$datatypeCMref${strain}${total}_empty --binSize 50 --sortUsing max --skipZeros -p=max --maxThreshold 100000
            	plotHeatmap -m $COV/$datatypeCMref${strain}${total}_empty -out $COV/${datatype}CMref${strain}${total}_emptyheatmap --colorList white,blue --sortUsing max
            	plotProfile -m $COV/$datatypeCMref${strain}${total}_empty -out $COV/${datatype}CMref${strain}${total}_emptyprofile --perGroup --plotType=fill
 			else
 			fi
 			
elif [ $datatype = "H3K27me3" or "H3K36me3" or "H3K56ac" ];
then
    for size in 250 500
    do
        total=$(expr $size \* 2)
        for strain in HEG4 EG4 A123 A119
        do
            BEDFILE=$TMP
            NEWBASE=${strain}_shared.${total}
            # I think this currently runs but it's ridiculously slow for some reason
            bamCoverage -b $BAM/$datatype/${strain}.rep1.bam -o $COV/$datatype${strain}.rep1.bw -bs=1 -p=max
            bamCoverage -b $BAM/$datatype/${strain}.rep2.bam -o $COV/$datatype${strain}.rep2.bw -bs=1 -p=max
         
            if [ $tmpfoldername = "04_background_sites" ]; then
           		
            	# Background_rep1
            	computeMatrix reference-point -S $COV/$datatype${strain}.rep1.bw -R tmp/04_background_sites/Simulate0001_${strain}_background_sites.bed   --upstream 2000 --downstream 2000 -out $COV/$datatypeCMref${strain}${total}rep1_background --binSize 50 --sortUsing max --skipZeros -p=max --maxThreshold 100000
            	plotHeatmap -m $COV/$datatypeCMref${strain}${total}rep1_background -out $COV/${datatype}CMref${strain}${total}rep1_backgroundheatmap --colorList white,blue --sortUsing max
            	plotProfile -m $COV/$datatypeCMref${strain}${total}rep1_background -out $COV/${datatype}CMref${strain}${total}rep1_backgroundprofile --perGroup --plotType=fill
           	 	# Background_rep2
            	computeMatrix reference-point -S $COV/$datatype${strain}.rep2.bw -R tmp/04_background_sites/Simulate0001_${strain}_background_sites.bed   --upstream 2000 --downstream 2000 -out $COV/$datatypeCMref${strain}${total}rep2_background --binSize 50 --sortUsing max --skipZeros -p=max --maxThreshold 100000
            	plotHeatmap -m $COV/$datatypeCMref${strain}${total}rep2_background -out $COV/${datatype}CMref${strain}${total}rep2_backgroundheatmap --colorList white,blue --sortUsing max
            	plotProfile -m $COV/$datatypeCMref${strain}${total}rep2_background -out $COV/${datatype}CMref${strain}${total}rep2_backgroundprofile --perGroup --plotType=fill    

			elif [ $tmpfoldername = "05_processed_shared"  ];

           		# Shared_rep1
            	computeMatrix reference-point -S $COV/$datatype${strain}.rep1.bw -R tmp/05_processed_shared/${strain}_shared_sites.bed   --upstream 2000 --downstream 2000 -out $COV/$datatypeCMref${strain}${total}rep1_shared --binSize 50 --sortUsing max --skipZeros -p=max --maxThreshold 100000
            	plotHeatmap -m $COV/$datatypeCMref${strain}${total}rep1_shared -out $COV/${datatype}CMref${strain}${total}rep1_sharedheatmap --colorList white,blue --sortUsing max
            	plotProfile -m $COV/$datatypeCMref${strain}${total}rep1_shared -out $COV/${datatype}CMref${strain}${total}rep1_sharedprofile --perGroup --plotType=fill
            	# Shared_rep2
            	computeMatrix reference-point -S $COV/$datatype${strain}.rep2.bw -R tmp/05_processed_shared/${strain}_shared_sites.bed   --upstream 2000 --downstream 2000 -out $COV/$datatypeCMref${strain}${total}rep2_shared --binSize 50 --sortUsing max --skipZeros -p=max --maxThreshold 100000
            	plotHeatmap -m $COV/$datatypeCMref${strain}${total}rep2_shared -out $COV/${datatype}CMref${strain}${total}rep2_sharedheatmap --colorList white,blue --sortUsing max
            	plotProfile -m $COV/$datatypeCMref${strain}${total}rep2_shared -out $COV/${datatype}CMref${strain}${total}rep2_sharedprofile --perGroup --plotType=fill
            
            
            elif [ $tmpfoldername = "03_processed_empty" ];

            	# Empty_rep1
            	computeMatrix reference-point -S $COV/$datatype${strain}.rep1.bw -R tmp/03_processed_empty/${strain}_empty_sites.bed    --upstream 2000 --downstream 2000 -out $COV/$datatypeCMref${strain}${total}rep1_empty --binSize 50 --sortUsing max --skipZeros -p=max --maxThreshold 100000
           		plotHeatmap -m $COV/$datatypeCMref${strain}${total}rep1_empty -out $COV/${datatype}CMref${strain}${total}rep1_emptyheatmap --colorList white,blue --sortUsing max
            	plotProfile -m $COV/$datatypeCMref${strain}${total}rep1_empty -out $COV/${datatype}CMref${strain}${total}rep1_emptyprofile --perGroup --plotType=fill
            	# Empty_rep2
            	computeMatrix reference-point -S $COV/$datatype${strain}.rep2.bw -R tmp/03_processed_empty/${strain}_empty_sites.bed    --upstream 2000 --downstream 2000 -out $COV/$datatypeCMref${strain}${total}rep2_empty --binSize 50 --sortUsing max --skipZeros -p=max --maxThreshold 100000
            	plotHeatmap -m $COV/$datatypeCMref${strain}${total}rep2_empty -out $COV/${datatype}CMref${strain}${total}rep2_emptyheatmap --colorList white,blue --sortUsing max
            	plotProfile -m $COV/$datatypeCMref${strain}${total}rep2_empty -out $COV/${datatype}CMref${strain}${total}rep2_emptyprofile --perGroup --plotType=fill

 			else
 			fi
else
fi
	done
done