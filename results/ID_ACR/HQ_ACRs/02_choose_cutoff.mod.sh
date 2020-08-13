#!/bin/bash -e

module load samtools
module load bedtools/2.28.0
module load ncbi-blast/2.9.0+

if [ $# -lt 5 ]; then
    echo "usage:sh cutoff.sh <density_file> <cutoff> <merge-gaps> <peaks-length> <reps_counts> <name> "
    echo "
    <density_file>: the output file from calculate_density step; 
    <cutoff>: the density cutoff to choose for each bin;
    <merge-gaps>: the maxium length allowed to merge nearby bins;
    <peaks-length>: the minium length of a peak;
    <reps_counts>: sample replication number;
    <genome_fasta>: PATH to the genome fasta file;
    <NCBI_MtPt_Db>: PATH to the NCBI blast+ database of plant mitochrodria and plastids sequences;
    <name>: the name of the output peak files
    the NCBI blast+ and BEDTools should be in the system path."
    exit -1; 
fi

density=$1
cutoff=$2
gap=$3
len=$4
reps=$5
name=$6
genome=Oryza_sativa.IRGSP-1.0.dna.toplevel_reformat.fa
path_MtPt_Db=/rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/ncbi-blast-2.8.1/Oryzasativa_Mitochloro.fa


CPUS=$SLURM_CPUS_ON_NODE
if [ ! $CPUS ]; then
    CPUS=1
fi
if [ ! -f  Oryzasativa_Mitochloro.fa.nhr ]; then
  makeblastdb -in Oryzasativa_Mitochloro.fa -dbtype nucl
fi

awk '$NF>('$reps'*'$cutoff')' $density |bedtools sort -i - |bedtools merge -d $gap -i - |awk '$3-$2>'$len'' >  $density.$cutoff.$gap.bedo

bedtools getfasta -fi $genome -bed $density.$cutoff.$gap.bedo -fo $density.$gap.fa
#export PATH=/rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/ncbi-blast-2.8.1
module load bedtools/2.28.0
module load ncbi-blast/2.9.0+

blastn -query $density.$gap.fa -out $name.mtpt -db /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/ncbi-blast-2.8.1/Oryzasativa_Mitochloro.fa -outfmt 7

sed '/#/d' $name.mtpt | cut -f1 |sed "s/\:/\t/g" |sed "s/-/\t/g" |bedtools sort -i - | uniq > $name.black

bedtools intersect -a  $density.$cutoff.$gap.bedo -b $name.black -v |cut -f1-3 > $name

rm $density.$gap.fa
rm $name.mtpt $density.$cutoff.$gap.bedo

