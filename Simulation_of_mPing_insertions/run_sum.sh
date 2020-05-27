#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --time=4:00:00
#SBATCH --output=run_sum.sh.%A_%a.stdout
#SBATCH -p intel
#SBATCH --chdir=./


#for head in Control Strains Shared Het Hom;
#for head in Ken09_shared Ken09_unique;
for head in Random_nofrq1 Random_nofrq2 Random_nofrq3
#for head in Random_frq
do
    echo "$head"
    /rhome/cjinfeng/BigData/software/bedtools2-2.19.0/bin/bedtools intersect -a $head\.gff -b MSU_r7.all.final.full.utr.gff3 -wao > $head\.intersect
    /rhome/cjinfeng/BigData/software/bedtools2-2.19.0/bin/bedtools closest -a $head\.gff -b MSU_r7.all.final.mRNA.gff -d > $head\.mRNA.intersect
    python mPing_position.py --input $head\.intersect --mrna $head\.mRNA.intersect
    python mPing_intron.py --input $head\.intersect
    python mPing_intergenic.py --input $head\.mRNA.intersect
    python mPing_position_2kb.py --input $head\.intersect --mrna $head\.mRNA.intersect
done

echo "Done"
