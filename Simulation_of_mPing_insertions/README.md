
## Generate frequency matrix from known mPing insertion sites


```
python pictogram_matrix.py ./pictogram/RIL.tsd.fa > ./pictogram/RIL.tsd.matrix

```

## Simulate random mPing insertions in rice genome

```
echo "Simulate 100 random mPings in genome without using frequency matrix (only use the matrix to determine strad)"
python mPing_sim_9merTSDV2_Random.py --input pictogram/somatic.tsd.matrix --output simulateV2_Random_TSD9mer_somaticMat --replicate 1 --size 100 --use_freq 0

echo "Simulate 100 random mPings in genome using frequency matrix (this will be slower than the above one)"
<<<<<<< HEAD
python mPing_sim_9merTSDV2_Random.py --input pictogram/somatic.tsd.matrix --output simulateV2_Random_TSD9mer_somaticMat --replicate 1 --size 100 --use_freq 1
=======
python mPing_sim_9merTSDV2_Random.py --input pictogram/somatic.tsd.matrix --output simulateV2_Random_TSD9mer_somaticMat --replicate 2 --size 100 --use_freq 1

echo "Or run this script"
sbatch mPing_sim_9merTSDV2_Random.sh
```


## Summary mPing distribution

```
echo "Simulated mPings without using frequency matrix"
cp simulateV2_Random_TSD9mer_somaticMat/Simulate0001.gff Random_nofrq.gff

echo "Simulated mPings using frequency matrix"
cp simulateV2_Random_TSD9mer_somaticMat/Simulate0002.gff Random_frq.gff

echo "summary"


>>>>>>> 29c91d6ec7eb68e0720155acc1e1388374169381
```
