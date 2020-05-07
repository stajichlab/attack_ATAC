#!/usr/bin/bash
#SBATCH -p short

module load bedtools
OUTDIR=windows
GENOMETAB=genome/Oryza_sativa.IRGSP-1.0.30.chroms.tab
if [ ! -f $GENOMETAB ]; then
	module load samtools
	samtools faidx genome/Oryza_sativa.IRGSP-1.0.30.dna.toplevel.fa
	cut -f1,2 genome/Oryza_sativa.IRGSP-1.0.30.dna.toplevel.fa.fai > $GENOMETAB
fi
mkdir -p $OUTDIR



bedtools makewindows -g genome/Oryza_sativa.IRGSP-1.0.30.chroms.tab -w 5000  > $OUTDIR/Oryza_sativa.IRGSP-1.0.30.windows_5kb.bed
bedtools makewindows -g genome/Oryza_sativa.IRGSP-1.0.30.chroms.tab -w 10000  > $OUTDIR/Oryza_sativa.IRGSP-1.0.30.windows_10kb.bed
bedtools makewindows -g genome/Oryza_sativa.IRGSP-1.0.30.chroms.tab -w 1000  > $OUTDIR/Oryza_sativa.IRGSP-1.0.30.windows_1kb.bed
