#BSUB -J benchmark
#BSUB -W 72:00
#BSUB -o  /work3/felteu/logs/graphpart/
#BSUB -e /work3/felteu/logs/graphpart
#BSUB -q hpc
#BSUB -n 12
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=40GB]"

for LENGTH in 50 100 150 200 250 300 350 400 450 500 550 600 650 700 850 900 950 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000
do
graphpart mmseqs2 -ff "runtime_benchmark/${LENGTH}_aas.fasta" \
--threshold 0.30 \
--labels-name label \
--out-file "runtime_benchmark/mmseqs2_${LENGTH}_aas_result.csv"

done


rm -rf runtime_benchmark/*.csv*
