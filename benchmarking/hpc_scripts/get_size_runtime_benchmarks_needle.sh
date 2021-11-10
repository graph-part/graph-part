#BSUB -J benchmark_rna
#BSUB -W 72:00
#BSUB -o  /work3/felteu/logs/graphpart/
#BSUB -e /work3/felteu/logs/graphpart
#BSUB -q hpc
#BSUB -n 12
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=40GB]"

CONDA_BASE=$(conda info --base) ; source $CONDA_BASE/etc/profile.d/conda.sh
conda activate gp-env

for NUMBER in 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000 13000 14000 15000 16000 17000 18000 19000 20000 21000 22000 23000 24000 25000 26000 27000 28000 29000 30000
do
graphpart needle -ff "/zhome/1d/8/153438/experiments/graph-part/benchmarking/runtime_benchmark/${NUMBER}_seqs.fasta" \
--threshold 0.30 \
--threads 24 \
--chunks 300 \
--labels-name label \
--triangular \
--out-file "/zhome/1d/8/153438/experiments/graph-part/benchmarking/runtime_benchmark/needle_${NUMBER}_seqs_result.csv"

done

rm -rf /zhome/1d/8/153438/experiments/graph-part/benchmarking/runtime_benchmark/*.csv*
