#!/usr/bin/bash
#SBATCH -p short --mem 8gb -N 1 -n 8 --out logs/background_sites.log

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
BACKGROUND=data/background
GENOME=genome/Oryza_sativa.IRGSP-1.0.30.dna.toplevel.fa
if [ ! -f $GENOME.fai ]; then
	module load samtools
	samtools faidx $GENOME
fi
# genome lengths file
cut -f 1,2 $GENOME.fai > $GENOME.genome

COV=coverage/processed_background
TMP=tmp/04_background_sites # I like to keep a clean folder so I would use tmp/05_processed_background
 # then you can designate in your .gitignore file tmp/*
 # and also have a cleanup step that rm -rf tmp

mkdir -p $COV $TMP
for infile in $(ls $BACKGROUND/*.gff)
do
    base=$(basename $infile .gff)
    background=$TMP/${base}_background_sites
    cut -f1,4,5,9 $infile | perl -p -e 's/^Chr//; s/ID=//;' > $background.bed

	for size in 250 500
	do
		total=$(expr $size \* 2)
		# use -b to add to both sides
		BEDFILE=$background.${total}_window.bed
		bedtools slop -i $background.bed -g $GENOME.genome -b $size  > $BEDFILE
		NEWBASE=$(basename $BEDFILE _window.bed)
		for strain in HEG4 EG4 A123 A119
		do
			parallel --rpl '{///} $Global::use{"File::Basename"} ||= eval "use File::Basename; 1;"; $_ = basename(dirname($_));' \
				-j $CPUS mkdir -p ${COV}/{///} \&\& mosdepth -t $MOSCPU -f $GENOME -x -n --by $BEDFILE $COV/{///}/{///}_{/.}.${NEWBASE} {} ::: $(ls $BAM/*/${strain}*.bam)
		done
	done
done
