#!/bin/bash
#SBATCH -p short

module load bcftools

BEDFILE="all_empty_sites.500_window.bed"

cat $BEDFILE | while read CHROM START END DETAILS
do
 #echo "Chrom=$CHROM START=$START END=$END"
 Chrom=$CHROM 
 START=$START 
 END=$END
 DETAILS=$DETAILS
 tabix allc_A119.tsv.gz Chr${Chrom}:${START}-${END} > "./newbedfiles/A119_Chr${Chrom}:${START}-${END}.bed"
 tabix allc_A123.tsv.gz Chr${Chrom}:${START}-${END} > "./newbedfiles/A123_Chr${Chrom}:${START}-${END}.bed"
 tabix allc_EGH.tsv.gz Chr${Chrom}:${START}-${END} > "./newbedfiles/EGH_Chr${Chrom}:${START}-${END}.bed"
 tabix allc_HEG4.tsv.gz Chr${Chrom}:${START}-${END} > "./newbedfiles/HEG4_Chr${Chrom}:${START}-${END}.bed"
 #tabix allc_A119.tsv.gz Chr1:9000-9100 > "./bedfiles/${Chrom}:${START}-${END}.bed"
done