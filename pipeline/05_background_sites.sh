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
TMP=05_processed_background
mkdir -p $COV $TMP

#
cat $MPING/Simulate0001.gff | cut -f 1,4,5,9 | > $TMP/Simulate0001_background_sites.bed
cat $MPING/Simulate0002.gff | cut -f 1,4,5,9 | > $TMP/Simulate0002_background_sites.bed


# use -b to add to both sides
bedtools slop -i $TMP/Simulate0001_background_sites.bed -g $GENOME.genome -b 250  > $TMP/Simulate0001_background_sites.500_window.bed
bedtools slop -i $TMP/Simulate0001_background_sites.bed -g $GENOME.genome -b 500  > $TMP/Simulate0001_background_sites.1000_window.bed

bedtools slop -i $TMP/Simulate0002_background_sites.bed -g $GENOME.genome -b 250  > $TMP/Simulate0002_background_sites.500_window.bed
bedtools slop -i $TMP/Simulate0002_background_sites.bed -g $GENOME.genome -b 500  > $TMP/Simulate0002_background_sites.1000_window.bed

for BED in $(ls $TMP/*_window.bed)
do
	NEWBASE=$(basename $BED .bed)
	for BAMFILE  in $(ls $BAM/*/*.bam)
	do
		EPI=ATAC
		mosdepth -t $CPUS -f $GENOME -n -x --by $BED $COV/${EPI}.${NEWBASE} $BAMFILE
	done
done

for BED in $(ls $TMP/*_window.bed)
do
	NEWBASE=$(basename $BED .bed)
	STRAIN=$(echo $NEWBASE | perl -p -e 's/_empty\.\S+//')
	# could replace this with parallel code
	parallel --rpl '{///} $Global::use{"File::Basename"} ||= eval "use File::Basename; 1;"; $_ = basename(dirname($_));' \
		 -j $CPUS mosdepth -f $GENOME -x -n --by $BED $COV/{///}_{/.}${NEWBASE} {} ::: $(ls $BAM/*/*.bam)
	#for BAMFILE  in $(ls $BAM/*/$STRAIN.*.bam)
	#do
	#	EPI=$(basename `dirname $BAMFILE`)
	#	SAMPLE=$(basename $BAMFILE .bam)
	#	mosdepth -t $CPUS -f $GENOME -n -x --by $BED $COV/${EPI}.${SAMPLE}.${NEWBASE} $BAMFILE
	#done
done

mkdir coverage/processed_background/tmp
mv coverage/processed_background/*bed.gz.csi coverage/processed_background/tmp
mv coverage/processed_background/*global.dist.txt coverage/processed_background/tmp
mv coverage/processed_background/*region.dist.txt coverage/processed_background/tmp
mv coverage/processed_background/*.summary.txt coverage/processed_background/tmp

mkdir coverage/processed_background/ATAC
mkdir coverage/processed_background/H3K27me3
mkdir coverage/processed_background/H3K36me3
mkdir coverage/processed_background/H3K56ac

mv coverage/processed_background/ATAC* coverage/processed_background/ATAC
mv coverage/processed_background/H3K27me3* coverage/processed_background/H3K27me3
mv coverage/processed_background/H3K36me3* coverage/processed_background/H3K36me3
mv coverage/processed_background/H3K56ac* coverage/processed_background/H3K56ac