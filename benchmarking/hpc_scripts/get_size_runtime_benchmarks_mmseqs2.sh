#BSUB -J benchmark
#BSUB -W 72:00
#BSUB -o  /work3/felteu/logs/graphpart/
#BSUB -e /work3/felteu/logs/graphpart
#BSUB -q hpc
#BSUB -n 12
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=40GB]"

for NUMBER in 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000 13000 14000 15000 16000 17000 18000 19000 20000 21000 22000 23000 24000 25000 26000 27000 28000 29000 30000 35000 40000 45000 50000 55000 65000 70000 75000 80000 85000 90000 95000 100000
do
graphpart mmseqs2 -ff "runtime_benchmark/${NUMBER}_seqs.fasta" \
--threshold 0.30 \
--labels-name label \
--out-file "runtime_benchmark/mmseqs2_${NUMBER}_seqs_result.csv"

done

rm -rf runtime_benchmark/*.csv*