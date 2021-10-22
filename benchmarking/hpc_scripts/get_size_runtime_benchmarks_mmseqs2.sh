#BSUB -J benchmark_rna
#BSUB -W 72:00
#BSUB -o  /work3/felteu/logs/graphpart/
#BSUB -e /work3/felteu/logs/graphpart
#BSUB -q hpc
#BSUB -n 12
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=40GB]"


conda activate gp-env

for NUMBER in 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000 13000 14000 15000 16000 17000 18000 19000 20000
do
graphpart mmseqs2 -ff "/zhome/1d/8/153438/experiments/graph-part/benchmarking/runtime_benchmark/${NUMBER}_seqs.fasta" \
--threshold 0.30 \
--labels-name label \
--out-file "/zhome/1d/8/153438/experiments/graph-part/benchmarking/runtime_benchmark/mmseqs2_${NUMBER}_seqs_result.csv"

done