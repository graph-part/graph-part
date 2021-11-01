mkdir results

##
## 1. Protein data
## 

# CD-HIT
python3 baselines/homology_partition_cd_hit.py --fasta-file data/protein/netgpi_dataset.fasta --threshold 0.3 --labels-name label --out-file results/protein_partition_cd_hit.csv
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/protein/netgpi_dataset.fasta --threshold 0.3 --labels-name label --out-file results/protein_reduction_cd_hit.csv

python3 partitioning_quality/align_partitions.py --partition-file results/protein_partition_cd_hit.csv --fasta-file data/protein/netgpi_dataset.fasta --out-file results/protein_partition_cd_hit_qc.csv
python3 partitioning_quality/align_partitions.py --partition-file results/protein_reduction_cd_hit.csv --fasta-file data/protein/netgpi_dataset.fasta --out-file results/protein_reduction_cd_hit_qc.csv

# MMseqs2
python3 baselines/homology_partition_mmseqs2.py --fasta-file data/protein/netgpi_dataset.fasta --threshold 0.3 --labels-name label --out-file results/protein_partition_mmseqs2_cmode0.csv --cluster-mode 0
python3 baselines/homology_partition_mmseqs2.py --fasta-file data/protein/netgpi_dataset.fasta --threshold 0.3 --labels-name label --out-file results/protein_partition_mmseqs2_cmode1.csv --cluster-mode 1
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/protein/netgpi_dataset.fasta --threshold 0.3 --labels-name label --out-file results/protein_reduction_mmseqs2_cmode0.csv --cluster-mode 0
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/protein/netgpi_dataset.fasta --threshold 0.3 --labels-name label --out-file results/protein_reduction_mmseqs2_cmode1.csv --cluster-mode 1

python3 partitioning_quality/align_partitions.py --partition-file results/protein_partition_mmseqs2_cmode0.csv --fasta-file data/protein/netgpi_dataset.fasta --out-file results/protein_partition_mmseqs2_cmode0_qc.csv
python3 partitioning_quality/align_partitions.py --partition-file results/protein_partition_mmseqs2_cmode1.csv --fasta-file data/protein/netgpi_dataset.fasta --out-file results/protein_partition_mmseqs2_cmode1_qc.csv
python3 partitioning_quality/align_partitions.py --partition-file results/protein_reduction_mmseqs2_cmode0.csv --fasta-file data/protein/netgpi_dataset.fasta --out-file results/protein_reduction_mmseqs2_cmode0_qc.csv
python3 partitioning_quality/align_partitions.py --partition-file results/protein_reduction_mmseqs2_cmode1.csv --fasta-file data/protein/netgpi_dataset.fasta --out-file results/protein_reduction_mmseqs2_cmode_1qc.csv

# Graph-Part
graphpart needle -ff data/protein/netgpi_dataset.fasta -th 0.3 -ln label -of results/protein_partition_graphpart_needle.csv -pa 5 --threads 12
graphpart mmseqs2 -ff data/protein/netgpi_dataset.fasta -th 0.3 -ln label -of results/protein_partition_graphpart_mmseqs2.csv -pa 5

python3 partitioning_quality/align_partitions.py --partition-file results/protein_partition_graphpart_needle.csv --fasta-file data/protein/netgpi_dataset.fasta --out-file results/protein_partition_graphpart_needle_qc.csv
python3 partitioning_quality/align_partitions.py --partition-file results/protein_partition_graphpart_mmseqs2.csv --fasta-file data/protein/netgpi_dataset.fasta --out-file results/protein_partition_graphpart_mmseqs2_qc.csv

##
## 2. RNA data
## 

# CD-HIT
python3 baselines/homology_partition_cd_hit.py --fasta-file data/rna/nRC_all.fasta --threshold 0.8 --labels-name label --out-file results/rna_partition_cd_hit.csv --nucleotide
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/rna/nRC_all.fasta --threshold 0.8 --labels-name label --out-file results/rna_reduction_cd_hit.csv --nucleotide

python3 partitioning_quality/align_partitions.py --partition-file results/rna_partition_cd_hit.csv --fasta-file data/rna/nRC_all.fasta --out-file results/rna_partition_cd_hit_qc.csv --nucleotide --matrix EDNAFULL --threads 20
python3 partitioning_quality/align_partitions.py --partition-file results/rna_reduction_cd_hit.csv --fasta-file data/rna/nRC_all.fasta --out-file results/rna_reduction_cd_hit_qc.csv --nucleotide --matrix EDNAFULL --threads 20

# MMseqs2
python3 baselines/homology_partition_mmseqs2.py --fasta-file data/rna/nRC_all.fasta --threshold 0.8 --labels-name label --out-file results/rna_partition_mmseqs2_cmode0.csv --cluster-mode 0
#python3 baselines/homology_partition_mmseqs2.py --fasta-file data/rna/nRC_all.fasta --threshold 0.8 --labels-name label --out-file results/rna_partition_mmseqs2_cmode1.csv --cluster-mode 1
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/rna/nRC_all.fasta --threshold 0.8 --labels-name label --out-file results/rna_reduction_mmseqs2_cmode0.csv --cluster-mode 0
#python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/rna/nRC_all.fasta --threshold 0.8 --labels-name label --out-file results/rna_reduction_mmseqs2_cmode1.csv --cluster-mode 1

python3 partitioning_quality/align_partitions.py --partition-file results/rna_partition_mmseqs2_cmode0.csv --fasta-file data/rna/nRC_all.fasta --out-file results/rna_partition_mmseqs2_cmode0_qc.csv --nucleotide --matrix EDNAFULL --threads 15
#python3 partitioning_quality/align_partitions.py --partition-file results/rna_partition_mmseqs2_cmode1.csv --fasta-file data/rna/nRC_all.fasta --out-file results/rna_partition_mmseqs2_cmode1_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --partition-file results/rna_reduction_mmseqs2_cmode0.csv --fasta-file data/rna/nRC_all.fasta --out-file results/rna_reduction_mmseqs2_cmode0_qc.csv --nucleotide --matrix EDNAFULL
#python3 partitioning_quality/align_partitions.py --partition-file results/rna_reduction_mmseqs2_cmode1.csv --fasta-file data/rna/nRC_all.fasta --out-file results/rna_reduction_mmseqs2_cmode1_qc.csv --nucleotide --matrix EDNAFULL

# Graph-Part
graphpart needle -ff data/rna/nRC_all.fasta -th 0.8 -ln label -of results/rna_partition_graphpart_needle.csv -pa 5 --threads 12 --nucleotide --matrix EDNAFULL
graphpart mmseqs2 -ff data/rna/nRC_all.fasta -th 0.8 -ln label -of results/rna_partition_graphpart_mmseqs2.csv -pa 5 --nucleotide

python3 partitioning_quality/align_partitions.py --partition-file results/rna_partition_graphpart_needle.csv --fasta-file data/rna/nRC_all.fasta --out-file results/rna_partition_graphpart_needle_qc.csv --nucleotide --matrix EDNAFULL --threads 20
python3 partitioning_quality/align_partitions.py --partition-file results/rna_partition_graphpart_mmseqs2.csv --fasta-file data/rna/nRC_all.fasta --out-file results/rna_partition_graphpart_mmseqs2_qc.csv --nucleotide --matrix EDNAFULL




graphpart needle -ff data/protein/DeepLoc/deeploc_all.fasta -th 0.3 -ln label -of results/deeploc_partition_graphpart_needle.csv -pa 5 --threads 20
graphpart needle -ff data/protein/SignalP5/train_set_reformatted.fasta -th 0.3 -ln label -of results/signalp_partition_graphpart_needle.csv -pa 5 --threads 20


##
## 3. other datasets
##


# DeepLoc
python3 baselines/homology_partition_cd_hit.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --threshold 0.3 --labels-name label --out-file results/deeploc_partition_cd_hit.csv
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --threshold 0.3 --labels-name label --out-file results/deeploc_reduction_cd_hit.csv

python3 baselines/homology_partition_mmseqs2.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --threshold 0.3 --labels-name label --out-file results/deeploc_partition_mmseqs2_cmode0.csv --cluster-mode 0
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --threshold 0.3 --labels-name label --out-file results/deeploc_reduction_mmseqs2_cmode0.csv --cluster-mode 0

# needle needs proper hpc job
graphpart mmseqs2 -ff data/protein/DeepLoc/deeploc_all.fasta -th 0.3 -ln label -of results/deeploc_partition_graphpart_mmseqs2.csv -pa 5

python3 partitioning_quality/align_partitions.py --partition-file results/deeploc_partition_cd_hit.csv --fasta-file data/protein/DeepLoc/deeploc_all.fasta --out-file results/deeploc_partition_cd_hit_qc.csv --nucleotide --matrix EDNAFULL --threads 20
python3 partitioning_quality/align_partitions.py --partition-file results/deeploc_reduction_cd_hit.csv --fasta-file data/protein/DeepLoc/deeploc_all.fasta --out-file results/deeploc_reduction_cd_hit_qc.csv --nucleotide --matrix EDNAFULL --threads 20

python3 partitioning_quality/align_partitions.py --partition-file results/deeploc_partition_mmseqs2_cmode0.csv --fasta-file data/protein/DeepLoc/deeploc_all.fasta --out-file results/deeploc_partition_mmseqs2_cmode0_qc.csv --matrix EDNAFULL --threads 15
python3 partitioning_quality/align_partitions.py --partition-file results/deeploc_reduction_mmseqs2_cmode0.csv --fasta-file data/protein/DeepLoc/deeploc_all.fasta --out-file results/deeploc_reduction_mmseqs2_cmode0_qc.csv --matrix EDNAFULL --threads 15

python3 partitioning_quality/align_partitions.py --partition-file results/deeploc_partition_graphpart_needle.csv --fasta-file data/protein/DeepLoc/deeploc_all.fasta --out-file results/deeploc_partition_graphpart_needle_qc.csv 
python3 partitioning_quality/align_partitions.py --partition-file results/deeploc_partition_graphpart_mmseqs2.csv --fasta-file data/protein/DeepLoc/deeploc_all.fasta --out-file results/deeploc_partition_graphpart_mmseqs2_qc.csv


# SignalP
python3 baselines/homology_partition_cd_hit.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --threshold 0.3 --labels-name label --out-file results/signalp_partition_cd_hit.csv
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --threshold 0.3 --labels-name label --out-file results/signalp_reduction_cd_hit.csv

python3 baselines/homology_partition_mmseqs2.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --threshold 0.3 --labels-name label --out-file results/signalp_partition_mmseqs2_cmode0.csv --cluster-mode 0
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --threshold 0.3 --labels-name label --out-file results/signalp_reduction_mmseqs2_cmode0.csv --cluster-mode 0

# needle needs proper hpc job
graphpart mmseqs2 -ff data/protein/SignalP/train_set_reformatted.fasta -th 0.3 -ln label -of results/signalp_partition_graphpart_mmseqs2.csv -pa 5

python3 partitioning_quality/align_partitions.py --partition-file results/signalp_partition_cd_hit.csv --fasta-file data/protein/SignalP/train_set_reformatted.fasta --out-file results/signalp_partition_cd_hit_qc.csv --threads 20
python3 partitioning_quality/align_partitions.py --partition-file results/signalp_reduction_cd_hit.csv --fasta-file data/protein/SignalP/train_set_reformatted.fasta --out-file results/signalp_reduction_cd_hit_qc.csv --threads 20

python3 partitioning_quality/align_partitions.py --partition-file results/signalp_partition_mmseqs2_cmode0.csv --fasta-file data/protein/SignalP/train_set_reformatted.fasta --out-file results/signalp_partition_mmseqs2_cmode0_qc.csv --threads 15
python3 partitioning_quality/align_partitions.py --partition-file results/signalp_reduction_mmseqs2_cmode0.csv --fasta-file data/protein/SignalP/train_set_reformatted.fasta --out-file results/signalp_reduction_mmseqs2_cmode0_qc.csv --threads 15

python3 partitioning_quality/align_partitions.py --partition-file results/signalp_partition_graphpart_needle.csv --fasta-file data/protein/SignalP/train_set_reformatted.fasta --out-file results/signalp_partition_graphpart_needle_qc.csv 
python3 partitioning_quality/align_partitions.py --partition-file results/signalp_partition_graphpart_mmseqs2.csv --fasta-file data/protein/SignalP/train_set_reformatted.fasta --out-file results/signalp_partition_graphpart_mmseqs2_qc.csv


# FLIP thermostability
python3 baselines/homology_partition_cd_hit.py --fasta-file data/protein/FLIP/raw_human_sequences.fasta --threshold 0.3 --labels-name label --out-file results/flip_human_partition_cd_hit.csv
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/protein/FLIP/raw_human_sequences.fasta --threshold 0.3 --labels-name label --out-file results/flip_human_reduction_cd_hit.csv

python3 baselines/homology_partition_mmseqs2.py --fasta-file data/protein/FLIP/raw_human_sequences.fasta --threshold 0.3 --labels-name label --out-file results/flip_human_partition_mmseqs2_cmode0.csv --cluster-mode 0
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/protein/FLIP/raw_human_sequences.fasta --threshold 0.3 --labels-name label --out-file results/flip_human_reduction_mmseqs2_cmode0.csv --cluster-mode 0

# needle needs proper hpc job
graphpart mmseqs2 -ff data/protein/FLIP/raw_human_sequences.fasta -th 0.3 -ln label -of results/flip_human_partition_graphpart_mmseqs2.csv -pa 5

python3 partitioning_quality/align_partitions.py --partition-file results/flip_human_partition_cd_hit.csv --fasta-file data/protein/FLIP/raw_human_sequences.fasta --out-file results/flip_human_partition_cd_hit_qc.csv --threads 64
python3 partitioning_quality/align_partitions.py --partition-file results/flip_human_reduction_cd_hit.csv --fasta-file data/protein/FLIP/raw_human_sequences.fasta --out-file results/flip_human_reduction_cd_hit_qc.csv --threads 64

python3 partitioning_quality/align_partitions.py --partition-file results/flip_human_partition_mmseqs2_cmode0.csv --fasta-file data/protein/FLIP/raw_human_sequences.fasta --out-file results/flip_human_partition_mmseqs2_cmode0_qc.csv --threads 64
python3 partitioning_quality/align_partitions.py --partition-file results/flip_human_reduction_mmseqs2_cmode0.csv --fasta-file data/protein/FLIP/raw_human_sequences.fasta --out-file results/flip_human_reduction_mmseqs2_cmode0_qc.csv  --threads 64

python3 partitioning_quality/align_partitions.py --partition-file results/flip_human_partition_graphpart_needle.csv --fasta-file data/protein/FLIP/raw_human_sequences.fasta --out-file results/flip_human_partition_graphpart_needle_qc.csv --threads 64
python3 partitioning_quality/align_partitions.py --partition-file results/flip_human_partition_graphpart_mmseqs2.csv --fasta-file data/protein/FLIP/raw_human_sequences.fasta --out-file results/flip_human_partition_graphpart_mmseqs2_qc.csv --threads 64






