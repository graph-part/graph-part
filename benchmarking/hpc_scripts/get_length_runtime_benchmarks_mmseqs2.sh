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

for LENGTH in 50 100 150 200 250 300 350 400 450 500 550 600 650 700 850 900 950 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000
do
graphpart mmseqs2 -ff "/zhome/1d/8/153438/experiments/graph-part/benchmarking/runtime_benchmark/${LENGTH}_aas.fasta" \
--threshold 0.30 \
--labels-name label \
--out-file "/zhome/1d/8/153438/experiments/graph-part/benchmarking/runtime_benchmark/mmseqs2_${LENGTH}_aas_result.csv"

done


rm -rf /zhome/1d/8/153438/experiments/graph-part/benchmarking/runtime_benchmark/*.csv*
