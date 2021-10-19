#BSUB -J benchmark_rna
#BSUB -W 72:00
#BSUB -o  /work3/felteu/logs/graphpart/
#BSUB -e /work3/felteu/logs/graphpart
#BSUB -q hpc
#BSUB -n 12
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=40GB]"


conda activate gp-env

graphpart -ff /zhome/1d/8/153438/experiments/graph-part/benchmarking/data/mloci_dataset.fasta \
--threshold 0.01 \
--threads 24 \
--chunks 300 \
-labels-name label \
--save-checkpoint-path /work3/felteu/rna_benchmark \
--nucleotide \
--matrix EDNAFULL\
--triangular