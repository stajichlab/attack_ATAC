#!/usr/bin/bash
#SBATCH -p short --mem 64gb -N 1 -n 24 --out logs/coverage.log

module load bedtools
module unload perl
module load parallel
module load mosdepth

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
GENOME=genome/Oryza_sativa.IRGSP-1.0.30.dna.toplevel.fa
COV=coverage/processed_windows
mkdir -p $COV
# could do this with parallel with 2 arguments but that's prob even more confusing...
for WINDOW in 10000 5000 1000 500
do
    parallel --rpl '{///} $Global::use{"File::Basename"} ||= eval "use File::Basename; 1;"; $_ = basename(dirname($_));' \
		 -j $CPUS mkdir -p ${COV}/{///} \&\& mosdepth -t $MOSCPU -f $GENOME -x -n --by $WINDOW $COV/{///}/{///}_{/.}_w${WINDOW} {} ::: $(ls $BAM/*/*.bam)
done
