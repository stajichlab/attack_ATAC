#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=4:00:00
#SBATCH --output=mPing_sim_9merTSDV2_Random.sh.%A_%a.stdout
#SBATCH -p intel
#SBATCH --chdir=./


start=`date +%s`

CPU=$SLURM_NTASKS
if [ ! $CPU ]; then
   CPU=2
fi

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=1
fi

echo "CPU: $CPU"
echo "N: $N"

echo "Simulate random mPings in genome without using frequency matrix (only use the matrix to determine strad)"
#python mPing_sim_9merTSDV2_Random.py --input pictogram/somatic.tsd.matrix --output simulateV2_Random_TSD9mer_somaticMat --replicate 1 --size 1000 --use_freq 0
#python mPing_sim_9merTSDV2_Random.py --input pictogram/somatic.tsd.matrix --output simulateV2_Random_TSD9mer_somaticMat --replicate 2 --size 1000 --use_freq 0
#python mPing_sim_9merTSDV2_Random.py --input pictogram/somatic.tsd.matrix --output simulateV2_Random_TSD9mer_somaticMat --replicate 3 --size 1000 --use_freq 0

echo "Simulate random mPings in genome using frequency matrix (this will be slower than the above one)"
#python mPing_sim_9merTSDV2_Random.py --input pictogram/somatic.tsd.matrix --output simulateV2_Random_TSD9mer_somaticMat --replicate 4 --size 1000 --use_freq 1
python mPing_sim_9merTSDV2_Random.py --input pictogram/somatic.tsd.matrix --output simulateV2_Random_TSD9mer_somaticMat --replicate 5 --size 1000 --use_freq 1
python mPing_sim_9merTSDV2_Random.py --input pictogram/somatic.tsd.matrix --output simulateV2_Random_TSD9mer_somaticMat --replicate 6 --size 1000 --use_freq 1

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
