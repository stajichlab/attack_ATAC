#!/usr/bin/bash
module load samtools
module load bedtools
mkdir -p genome logs
pushd genome
FILE=Oryza_sativa.IRGSP-1.0.30.dna.toplevel.fa
if [[ ! -s $FILE ]]; then
	curl -L ftp://ftp.ensemblgenomes.org/pub/release-30/plants/fasta/oryza_sativa/dna/$FILE.gz | pigz -d > $FILE
fi
# index for samtools retrieval
samtools faidx $FILE
# make a genome file for bedtools
cut -f1,2 $FILE.fai > Oryza_sativa.IRGSP-1.0.30.chroms.tab
awk '{OFS="\t"} {print $1,1,$2}' $FILE.fai > Oryza_sativa.IRGSP-1.0.30.chroms.bed

FILE=Oryza_sativa.IRGSP-1.0.30.gff3
if [[ ! -s $FILE ]]; then
	curl -L ftp://ftp.ensemblgenomes.org/pub/release-30/plants/gff3/oryza_sativa/$FILE.gz | pigz -d > $FILE
fi

# this only works on UCR cluster
rsync -a /bigdata/stajichlab/jstajich/projects/attack_ATAC/bam ./
