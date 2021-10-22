#BSUB -J benchmark_rna
#BSUB -W 72:00
#BSUB -o  /work3/felteu/logs/graphpart/
#BSUB -e /work3/felteu/logs/graphpart
#BSUB -q hpc
#BSUB -n 12
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=40GB]"


conda activate gp-env

for LENGTH in 50 100 150 200 250 300 350 400 450 500
do
graphpart mmseqs2 -ff "/zhome/1d/8/153438/experiments/graph-part/benchmarking/runtime_benchmark/${LENGTH}_aas.fasta" \
--threshold 0.30 \
--labels-name label \
--out-file "/zhome/1d/8/153438/experiments/graph-part/benchmarking/runtime_benchmark/mmseqs2_${LENGTH}_result.csv"

done

