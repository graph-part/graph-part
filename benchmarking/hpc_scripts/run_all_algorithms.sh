#!/bin/bash
#SBATCH --partition=compute
#SBATCH --ntasks=1 --cpus-per-task=64 --mem=40G
#SBATCH --time=84:00:00
#SBATCH --output=/novo/users/fegt/graph-part/benchmarking/logs/%j.out
#SBATCH --error=/novo/users/fegt/graph-part/benchmarking/logs/%j.err

lscpu | egrep 'Model name|Socket|Thread|NUMA|CPU\(s\)'

source /novo/users/fegt/.bashrc
conda activate graphpart
cd /novo/users/fegt/graph-part/benchmarking

pwd

# GraphPart

graphpart needle -ff data/protein/NetGPI/netgpi_dataset.fasta -th 0.3 -ln label -of results/netgpi_partition_graphpart_needle.csv -pa 5 --threads 64
graphpart mmseqs2 -ff data/protein/NetGPI/netgpi_dataset.fasta -th 0.3 -ln label -of results/netgpi_partition_graphpart_mmseqs2.csv -pa 5

graphpart needle -ff data/rna/nRC/nRC_all.fasta -th 0.8 -ln label -of results/nrc_partition_graphpart_needle.csv -pa 5 --threads 64 --nucleotide --matrix EDNAFULL
graphpart mmseqs2 -ff data/rna/nRC/nRC_all.fasta -th 0.8 -ln label -of results/nrc_partition_graphpart_mmseqs2.csv -pa 5 --nucleotide

graphpart needle -ff data/protein/DeepLoc/deeploc_all.fasta -th 0.3 -ln label -of results/deeploc_partition_graphpart_needle.csv -pa 5 --threads 54
graphpart mmseqs2 -ff data/protein/DeepLoc/deeploc_all.fasta -th 0.3 -ln label -of results/deeploc_partition_graphpart_mmseqs2.csv -pa 5

graphpart needle -ff data/protein/SignalP/train_set_reformatted.fasta -th 0.3 -ln label -of results/signalp_partition_graphpart_needle.csv -pa 5 --threads 64
graphpart mmseqs2 -ff data/protein/SignalP/train_set_reformatted.fasta -th 0.3 -ln label -of results/signalp_partition_graphpart_mmseqs2.csv -pa 5

graphpart needle -ff data/dna/DeePromoter/deepromoter_tata_human.fasta -th 0.8 -of results/deepromoter_partition_graphpart_needle.csv -pa 5 --threads 64 --nucleotide --matrix EDNAFULL
graphpart mmseqs2 -ff data/dna/DeePromoter/deepromoter_tata_human.fasta -th 0.8 -of results/deepromoter_partition_graphpart_mmseqs2.csv -pa 5 --nucleotide

graphpart needle -ff data/dna/Histone/H3.fasta -th 0.8 -of results/histone_partition_graphpart_needle.csv -pa 5 --threads 64 --nucleotide --matrix EDNAFULL --labels-name label
graphpart mmseqs2 -ff data/dna/Histone/H3.fasta -th 0.8 -of results/histone_partition_graphpart_mmseqs2.csv -pa 5 --nucleotide --labels-name label

graphpart needle -ff data/rna/SPOT-RNA-1D/spot_rna_1d.fasta -th 0.8 -of results/spotrna1d_partition_graphpart_needle.csv -pa 5 --threads 64 --nucleotide --matrix EDNAFULL --labels-name label
graphpart mmseqs2 -ff data/rna/SPOT-RNA-1D/spot_rna_1d.fasta -th 0.8 -of results/spotrna1d_partition_graphpart_mmseqs2.csv -pa 5 --nucleotide --labels-name label

graphpart needle -ff data/protein/NetSurfP/netsurfp.fasta -th 0.25 -ln label -of results/netsurfp_partition_graphpart_needle.csv -pa 5 --threads 64
graphpart mmseqs2 -ff data/protein/NetSurfP/netsurfp.fasta -th 0.25 -ln label -of results/netsurfp_partition_graphpart_mmseqs2.csv -pa 5

# CD-HIT partitioning

python3 baselines/homology_partition_cd_hit.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --threshold 0.3 --out-file results/netgpi_partition_cdhit.csv --labels-name label
python3 baselines/homology_partition_cd_hit.py --fasta-file data/rna/nRC/nRC_all.fasta --threshold 0.8 --out-file results/nrc_partition_cdhit.csv --nucleotide --labels-name label
python3 baselines/homology_partition_cd_hit.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --threshold 0.3 --out-file results/deeploc_partition_cdhit.csv --labels-name label
python3 baselines/homology_partition_cd_hit.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --threshold 0.3 --out-file results/signalp_partition_cdhit.csv --labels-name label
python3 baselines/homology_partition_cd_hit.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --threshold 0.8 --out-file results/deepromoter_partition_cdhit.csv --nucleotide
python3 baselines/homology_partition_cd_hit.py --fasta-file data/dna/Histone/H3.fasta --threshold 0.8 --out-file results/histone_partition_cdhit.csv --nucleotide --labels-name label
python3 baselines/homology_partition_cd_hit.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --threshold 0.8 --out-file results/spotrna1d_partition_cdhit.csv --nucleotide --labels-name label
python3 baselines/homology_partition_cd_hit.py --fasta-file data/protein/NetSurfP/netsurfp.fasta --threshold 0.25 --out-file results/netsurfp_partition_cdhit.csv

# MMseqs2 partitioning

python3 baselines/homology_partition_mmseqs2.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --threshold 0.3 --out-file results/netgpi_partition_mmseqs2.csv --labels-name label
python3 baselines/homology_partition_mmseqs2.py --fasta-file data/rna/nRC/nRC_all.fasta --threshold 0.8 --out-file results/nrc_partition_mmseqs2.csv --labels-name label
python3 baselines/homology_partition_mmseqs2.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --threshold 0.3 --out-file results/deeploc_partition_mmseqs2.csv --labels-name label
python3 baselines/homology_partition_mmseqs2.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --threshold 0.3 --out-file results/signalp_partition_mmseqs2.csv --labels-name label
python3 baselines/homology_partition_mmseqs2.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --threshold 0.8 --out-file results/deepromoter_partition_mmseqs2.csv
python3 baselines/homology_partition_mmseqs2.py --fasta-file data/dna/Histone/H3.fasta --threshold 0.8 --out-file results/histone_partition_mmseqs2.csv --labels-name label
python3 baselines/homology_partition_mmseqs2.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --threshold 0.8 --out-file results/spotrna1d_partition_mmseqs2.csv --labels-name label
python3 baselines/homology_partition_mmseqs2.py --fasta-file data/protein/NetSurfP/netsurfp.fasta --threshold 0.25 --out-file results/netsurfp_partition_mmseqs2.csv

# CD-HIT reduction

python3 baselines/homology_reduce_cd_hit.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --threshold 0.3 --out-file results/netgpi_reduction_cdhit.csv --labels-name label
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/rna/nRC/nRC_all.fasta --threshold 0.8 --out-file results/nrc_reduction_cdhit.csv --nucleotide --labels-name label
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --threshold 0.3 --out-file results/deeploc_reduction_cdhit.csv --labels-name label
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --threshold 0.3 --out-file results/signalp_reduction_cdhit.csv --labels-name label
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --threshold 0.8 --out-file results/deepromoter_reduction_cdhit.csv --nucleotide
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/dna/Histone/H3.fasta --threshold 0.8 --out-file results/histone_reduction_cdhit.csv --nucleotide --labels-name label
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --threshold 0.8 --out-file results/spotrna1d_reduction_cdhit.csv --nucleotide --labels-name label
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/protein/NetSurfP/netsurfp.fasta --threshold 0.25 --out-file results/netsurfp_reduction_cdhit.csv

# MMseqs2 reduction

python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --threshold 0.3 --out-file results/netgpi_reduction_mmseqs2.csv --labels-name label
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/rna/nRC/nRC_all.fasta --threshold 0.8 --out-file results/nrc_reduction_mmseqs2.csv --labels-name label
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --threshold 0.3 --out-file results/deeploc_reduction_mmseqs2.csv --labels-name label
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --threshold 0.3 --out-file results/signalp_reduction_mmseqs2.csv --labels-name label
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --threshold 0.8 --out-file results/deepromoter_reduction_mmseqs2.csv
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/dna/Histone/H3.fasta --threshold 0.8 --out-file results/histone_reduction_mmseqs2.csv --labels-name label
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --threshold 0.8 --out-file results/spotrna1d_reduction_mmseqs2.csv --labels-name label
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/protein/NetSurfP/netsurfp.fasta --threshold 0.25 --out-file results/netsurfp_reduction_mmseqs2.csv


# Sklearn random partitioning
python3 baselines/random_partition_sklearn.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --out-file results/netgpi_sklearn.csv --labels-name label
python3 baselines/random_partition_sklearn.py --fasta-file data/rna/nRC/nRC_all.fasta --out-file results/nrc_sklearn.csv --labels-name label
python3 baselines/random_partition_sklearn.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --out-file results/deeploc_sklearn.csv --labels-name label
python3 baselines/random_partition_sklearn.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --out-file results/signalp_sklearn.csv --labels-name label
python3 baselines/random_partition_sklearn.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --out-file results/deepromoter_sklearn.csv
python3 baselines/random_partition_sklearn.py --fasta-file data/dna/Histone/H3.fasta --out-file results/histone_sklearn.csv --labels-name label
python3 baselines/random_partition_sklearn.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --out-file results/spotrna1d_sklearn.csv --labels-name label
python3 baselines/random_partition_sklearn.py --fasta-file data/protein/NetSurfP/netsurfp.fasta --out-file results/netsurfp_sklearn.csv


# Hobohm 2 reduction
python3 baselines/homology_reduce_hobohm2.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --out-file results/netgpi_reduction_hobohm2.csv --labels-name label -al needle -nt 64
python3 baselines/homology_reduce_hobohm2.py --fasta-file data/rna/nRC/nRC_all.fasta --out-file results/nrc_reduction_hobohm2.csv --labels-name label -al needle -nt 64 --nucleotide -th 0.8
python3 baselines/homology_reduce_hobohm2.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --out-file results/deeploc_reduction_hobohm2.csv --labels-name label -al needle -nt 64
python3 baselines/homology_reduce_hobohm2.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --out-file results/signalp_reduction_hobohm2.csv --labels-name label -al needle -nt 64
python3 baselines/homology_reduce_hobohm2.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --out-file results/deepromoter_reduction_hobohm2.csv -al needle -nt 64 --nucleotide -th 0.8
python3 baselines/homology_reduce_hobohm2.py --fasta-file data/dna/Histone/H3.fasta --out-file results/histone_reduction_hobohm2.csv --labels-name label -al needle -nt 64 --nucleotide -th 0.8
python3 baselines/homology_reduce_hobohm2.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --out-file results/spotrna1d_reduction_hobohm2.csv --labels-name label -al needle -nt 64 --nucleotide -th 0.8
python3 baselines/homology_reduce_hobohm2.py --fasta-file data/protein/NetSurfP/netsurfp.fasta --out-file results/netsurfp_reduction_hobohm2.csv -al needle -nt 64 -th 0.25


python3 baselines/homology_reduce_hobohm2.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --out-file results/netgpi_reduction_hobohm2_mmseqs2.csv --labels-name label -al mmseqs2 -nt 64
python3 baselines/homology_reduce_hobohm2.py --fasta-file data/rna/nRC/nRC_all.fasta --out-file results/nrc_reduction_hobohm2_mmseqs2.csv --labels-name label -al mmseqs2 -nt 64 --nucleotide -th 0.8
python3 baselines/homology_reduce_hobohm2.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --out-file results/deeploc_reduction_hobohm2_mmseqs2.csv --labels-name label -al mmseqs2 -nt 64
python3 baselines/homology_reduce_hobohm2.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --out-file results/signalp_reduction_hobohm2_mmseqs2.csv --labels-name label -al mmseqs2 -nt 64
python3 baselines/homology_reduce_hobohm2.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --out-file results/deepromoter_reduction_hobohm2_mmseqs2.csv -al mmseqs2 -nt 64 --nucleotide -th 0.8
python3 baselines/homology_reduce_hobohm2.py --fasta-file data/dna/Histone/H3.fasta --out-file results/histone_reduction_hobohm2_mmseqs2.csv --labels-name label -al mmseqs2 -nt 64 --nucleotide -th 0.8
python3 baselines/homology_reduce_hobohm2.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --out-file results/spotrna1d_reduction_hobohm2_mmseqs2.csv --labels-name label -al mmseqs2 -nt 64 --nucleotide -th 0.8
python3 baselines/homology_reduce_hobohm2.py --fasta-file data/protein/NetSurfP/netsurfp.fasta --out-file results/netsurfp_reduction_hobohm2_mmseqs2.csv -al mmseqs2 -nt 64 -th 0.25