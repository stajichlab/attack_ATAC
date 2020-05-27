#!/usr/bin/bash
#SBATCH -p short --mem 8gb -N 1 -n 8 --out logs/coverage.log

module load bedtools
module load mosdepth
module unload perl
module load parallel


CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
	CPU=1
fi
CPUS=$CPU

BAM=bam
MPING=data/background
GENOME=genome/Oryza_sativa.IRGSP-1.0.30.dna.toplevel.fa
if [ ! -f $GENOME.fai ]; then
	module load samtools
	samtools faidx $GENOME
fi
# genome lengths file
cut -f 1,2 $GENOME.fai > $GENOME.genome

COV=coverage/processed_background
TMP=tmp/05_processed_background # I like to keep a clean folder so I would use tmp/05_processed_background
 # then you can designate in your .gitignore file tmp/*
 # and also have a cleanup step that rm -rf tmp

mkdir -p $COV $TMP
for infile in $(ls $MPING/*.gff)
do
    base=$(basename $infile .gff)
    background=$TMP/${base}_background_sites
    cut -f1,4,5,9 $infile | perl -p -e 's/^Chr//; s/ID=//;' > $background.bed

    # use -b to add to both sides
    bedtools slop -i $background.bed -g $GENOME.genome -b 250  > $background.1000_window.bed
    bedtools slop -i $background.bed -g $GENOME.genome -b 500  > $background.500_window.bed
done

for BED in $(ls $TMP/*_window.bed)
do
	NEWBASE=$(basename $BED .bed)
	for BAMFILE in $(ls $BAM/*/*.bam)
	do
		EPI=$(basename `dirname $BAMFILE`)
		mkdir -p $COV/$EPI
		mosdepth -t $CPUS -f $GENOME -n -x --by $BED $COV/${EPI}/${EPI}.${NEWBASE} $BAMFILE
	done
done
