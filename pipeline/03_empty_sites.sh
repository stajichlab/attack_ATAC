#!/usr/bin/bash
#SBATCH -p short --mem 8gb -N 1 -n 8 --out logs/coverage.log

module load bedtools
module load mosdepth

CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
	CPU=1
fi
CPUS=$CPU

BAM=bam
MPING=data/mPing
GENOME=genome/Oryza_sativa.IRGSP-1.0.30.dna.toplevel.fa
if [ ! -f $GENOME.fai ]; then
	module load samtools
	samtools faidx $GENOME
fi
# genome lengths file
cut -f 1,2 $GENOME.fai > $GENOME.genome

COV=coverage/pre_mPing
TMP=tmp
mkdir -p $COV $TMP
# Get HEG4 empty sites (mPing in EG4 but not in HEG4)
bedtools intersect -a $MPING/EG4.gff -b $MPING/HEG4.gff -v | cut -f 1,4,5,9 | perl -p -e 's/^Chr//; s/ID=//' > $TMP/HEG4_empty_sites.bed
# an vice versa
bedtools intersect -a $MPING/HEG4.gff -b $MPING/EG4.gff -v | cut -f 1,4,5,9 | perl -p -e 's/^Chr//; s/ID=//' > $TMP/EG4_empty_sites.bed

bedtools slop -i $TMP/EG4_empty_sites.bed -g $GENOME.genome -b 250  > $TMP/EG4_empty.500_window.bed
bedtools slop -i $TMP/EG4_empty_sites.bed -g $GENOME.genome -b 500  > $TMP/EG4_empty.1000_window.bed

bedtools slop -i $TMP/HEG4_empty_sites.bed -g $GENOME.genome -b 250 > $TMP/HEG4_empty.500_window.bed
bedtools slop -i $TMP/HEG4_empty_sites.bed -g $GENOME.genome -b 250 > $TMP/HEG4_empty.1000_window.bed

for BED in $(ls $TMP/*_window.bed)
do
	NEWBASE=$(basename $BED .bed)
	# could replace this with parallel code
	# parallel --rpl '{///} $Global::use{"File::Basename"} ||= eval "use File::Basename; 1;"; $_ = basename(dirname($_));' \
	#	 -j $CPUS mosdepth -f $GENOME -x -n --by $BED $COV/{///}_{/.}${NEWBASE} {} ::: $(ls $BAM/*/*.bam)
	for BAMFILE  in $(ls $BAM/*/*.bam)
	do
		EPI=$(basename `dirname $BAMFILE`)
		mosdepth -t $CPUS -f $GENOME -n -x --by $BED $COV/${EPI}.${NEWBASE} $BAMFILE
	done
done
