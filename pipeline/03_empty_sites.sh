#!/usr/bin/bash
#SBATCH -p short --mem 24gb -N 1 -n 24 --out logs/empty_sites.log

module load bedtools
module load mosdepth
module unload perl
module load parallel

CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
	CPU=1
fi

CPUS=$(expr $CPU / 2)
MOSCPU=2 # lets run mosdepth with 2 cores
if [[ $CPUS -lt 1 ]]; then
	CPUS=1
	MOSCPU=1
fi

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
NIPPO=tmp/03_processed_empty
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
for size in 250 500
do
	total=$(expr $size \* 2)
	for strain in HEG4 EG4 A123 A119
	do
			BEDFILE=$TMP/${strain}_empty.${total}_window.bed
			NEWBASE=${strain}_empty.${total}
			bedtools slop -i $TMP/${strain}_empty_sites.bed -g $GENOME.genome -b $size  > $BEDFILE

			parallel --rpl '{///} $Global::use{"File::Basename"} ||= eval "use File::Basename; 1;"; $_ = basename(dirname($_));' \
			-j $CPUS mkdir -p ${COV}/{///} \&\& mosdepth -t $MOSCPU -f $GENOME -x -n --by $BEDFILE $COV/{///}/{///}_{/.}.${NEWBASE} {} ::: $(ls $BAM/*/${strain}*.bam)

	done
done
