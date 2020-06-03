#!/usr/bin/bash
#SBATCH -p short --mem 24gb -N 1 -n 24 --out logs/shared_sites.log

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

COV=coverage/processed_shared
TMP=tmp/05_processed_shared
mkdir -p $COV $TMP

# Get shared sites
bedtools intersect -a $MPING/EG4.gff -b $MPING/HEG4.gff | cut -f 1,4,5,9 | perl -p -e 's/^Chr//; s/ID=//' > $TMP/HEG4_shared_sites.bed
bedtools intersect -a $MPING/EG4.gff -b $MPING/HEG4.gff | cut -f 1,4,5,9 | perl -p -e 's/^Chr//; s/ID=//' > $TMP/EG4_shared_sites.bed

bedtools intersect -a $MPING/A119.gff -b $MPING/A123.gff -v | cut -f 1,4,5,9 | perl -p -e 's/^Chr//; s/ID=//' > $TMP/A119_shared_sites.bed
bedtools intersect -a $MPING/A119.gff -b $MPING/A123.gff -v | cut -f 1,4,5,9 | perl -p -e 's/^Chr//; s/ID=//' > $TMP/A123_shared_sites.bed

# use -b to add to both sides
for size in 250 500
do
	total=$(expr $size \* 2)
	for strain in HEG4 EG4 A123 A119
	do
		BEDFILE=$TMP/${strain}_shared.${total}_window.bed
		NEWBASE=${strain}_shared.${total}
		bedtools slop -i $TMP/${strain}_shared_sites.bed -g $GENOME.genome -b $size  > $BEDFILE

		parallel --rpl '{///} $Global::use{"File::Basename"} ||= eval "use File::Basename; 1;"; $_ = basename(dirname($_));' \
		-j $CPUS mkdir -p ${COV}/{///} \&\& mosdepth -t $MOSCPU -f $GENOME -x -n --by $BEDFILE $COV/{///}/{///}_{/.}.${NEWBASE} {} ::: $(ls $BAM/*/${strain}*.bam)
	done
done