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

COV=coverage/processed_empty
TMP=tmp/03_processed_empty
mkdir -p $COV $TMP

# For HEG4 and EG4
# Get HEG4 empty sites (mPing in EG4 but not in HEG4)
bedtools intersect -a $MPING/EG4.gff -b $MPING/HEG4.gff -v | cut -f 1,4,5,9 | perl -p -e 's/^Chr//; s/ID=//' > $TMP/HEG4_empty_sites.bed
# an vice versa
bedtools intersect -a $MPING/HEG4.gff -b $MPING/EG4.gff -v | cut -f 1,4,5,9 | perl -p -e 's/^Chr//; s/ID=//' > $TMP/EG4_empty_sites.bed


# For A123 and A119
# Get A123 empty sites (mPing in A119 but not in A123)
bedtools intersect -a $MPING/A119.gff -b $MPING/A123.gff -v | cut -f 1,4,5,9 | perl -p -e 's/^Chr//; s/ID=//' > $TMP/A123_empty_sites.bed
# an vice versa
bedtools intersect -a $MPING/A123.gff -b $MPING/A119.gff -v | cut -f 1,4,5,9 | perl -p -e 's/^Chr//; s/ID=//' > $TMP/A119_empty_sites.bed

# use -b to add to both sides
bedtools slop -i $TMP/EG4_empty_sites.bed -g $GENOME.genome -b 250  > $TMP/EG4_empty.500_window.bed
bedtools slop -i $TMP/EG4_empty_sites.bed -g $GENOME.genome -b 500  > $TMP/EG4_empty.1000_window.bed

bedtools slop -i $TMP/HEG4_empty_sites.bed -g $GENOME.genome -b 250 > $TMP/HEG4_empty.500_window.bed
bedtools slop -i $TMP/HEG4_empty_sites.bed -g $GENOME.genome -b 500 > $TMP/HEG4_empty.1000_window.bed

bedtools slop -i $TMP/A123_empty_sites.bed -g $GENOME.genome -b 250  > $TMP/A123_empty.500_window.bed
bedtools slop -i $TMP/A123_empty_sites.bed -g $GENOME.genome -b 500  > $TMP/A123_empty.1000_window.bed

bedtools slop -i $TMP/A119_empty_sites.bed -g $GENOME.genome -b 250 > $TMP/A119_empty.500_window.bed
bedtools slop -i $TMP/A119_empty_sites.bed -g $GENOME.genome -b 500 > $TMP/A119_empty.1000_window.bed

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
	# parallel --rpl '{///} $Global::use{"File::Basename"} ||= eval "use File::Basename; 1;"; $_ = basename(dirname($_));' \
	#	 -j $CPUS mosdepth -f $GENOME -x -n --by $BED $COV/{///}_{/.}${NEWBASE} {} ::: $(ls $BAM/*/*.bam)
	for BAMFILE  in $(ls $BAM/*/$STRAIN.*.bam)
	do
		EPI=$(basename `dirname $BAMFILE`)
		SAMPLE=$(basename $BAMFILE .bam)
		mosdepth -t $CPUS -f $GENOME -n -x --by $BED $COV/${EPI}.${SAMPLE}.${NEWBASE} $BAMFILE
	done
done

mkdir coverage/processed_empty/tmp
mv coverage/processed_empty/*bed.gz.csi coverage/processed_empty/tmp
mv coverage/processed_empty/*global.dist.txt coverage/processed_empty/tmp
mv coverage/processed_empty/*region.dist.txt coverage/processed_empty/tmp

mkdir coverage/processed_empty/ATAC
mkdir coverage/processed_empty/H3K27me3
mkdir coverage/processed_empty/H3K36me3
mkdir coverage/processed_empty/H3K56ac

mv coverage/processed_empty/ATAC* coverage/processed_empty/ATAC
mv coverage/processed_empty/H3K27me3* coverage/processed_empty/H3K27me3
mv coverage/processed_empty/H3K36me3* coverage/processed_empty/H3K36me3
mv coverage/processed_empty/H3K56ac* coverage/processed_empty/H3K56ac
