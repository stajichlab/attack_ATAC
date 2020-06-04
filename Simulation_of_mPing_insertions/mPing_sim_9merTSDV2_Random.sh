#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=96G
#SBATCH --time=8:00:00
#SBATCH --output=mPing_sim_9merTSDV2_Random.sh.%A_%a.stdout
#SBATCH -p intel

module unload miniconda3
module load miniconda2

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
python mPing_sim_9merTSDV2_Random.py --input pictogram/somatic.tsd.matrix --output simulateV2_Random_TSD9mer_somaticMat --replicate $N --size 1000 --use_freq 0

echo "Simulate random mPings in genome using frequency matrix (this will be slower than the above one)"
python mPing_sim_9merTSDV2_Random.py --input pictogram/somatic.tsd.matrix --output simulateV2_Random_TSD9mer_somaticMat_wFreq --replicate $N --size 1000 --use_freq 1

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
