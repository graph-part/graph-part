mkdir results

##
## 1. NetGPI data
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
## 2. nRC data
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



##
## 3. DeepLoc data
##
python3 baselines/homology_partition_cd_hit.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --threshold 0.3 --labels-name label --out-file results/deeploc_partition_cd_hit.csv
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --threshold 0.3 --labels-name label --out-file results/deeploc_reduction_cd_hit.csv

python3 baselines/homology_partition_mmseqs2.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --threshold 0.3 --labels-name label --out-file results/deeploc_partition_mmseqs2_cmode0.csv --cluster-mode 0
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --threshold 0.3 --labels-name label --out-file results/deeploc_reduction_mmseqs2_cmode0.csv --cluster-mode 0

# needle needs proper hpc job
graphpart needle -ff data/protein/DeepLoc/deeploc_all.fasta -th 0.3 -ln label -of results/deeploc_partition_graphpart_needle.csv -pa 5 --threads 20
graphpart mmseqs2 -ff data/protein/DeepLoc/deeploc_all.fasta -th 0.3 -ln label -of results/deeploc_partition_graphpart_mmseqs2.csv -pa 5

python3 partitioning_quality/align_partitions.py --partition-file results/deeploc_partition_cd_hit.csv --fasta-file data/protein/DeepLoc/deeploc_all.fasta --out-file results/deeploc_partition_cd_hit_qc.csv --nucleotide --matrix EDNAFULL --threads 20
python3 partitioning_quality/align_partitions.py --partition-file results/deeploc_reduction_cd_hit.csv --fasta-file data/protein/DeepLoc/deeploc_all.fasta --out-file results/deeploc_reduction_cd_hit_qc.csv --nucleotide --matrix EDNAFULL --threads 20

python3 partitioning_quality/align_partitions.py --partition-file results/deeploc_partition_mmseqs2_cmode0.csv --fasta-file data/protein/DeepLoc/deeploc_all.fasta --out-file results/deeploc_partition_mmseqs2_cmode0_qc.csv --matrix EDNAFULL --threads 15
python3 partitioning_quality/align_partitions.py --partition-file results/deeploc_reduction_mmseqs2_cmode0.csv --fasta-file data/protein/DeepLoc/deeploc_all.fasta --out-file results/deeploc_reduction_mmseqs2_cmode0_qc.csv --matrix EDNAFULL --threads 15

python3 partitioning_quality/align_partitions.py --partition-file results/deeploc_partition_graphpart_needle.csv --fasta-file data/protein/DeepLoc/deeploc_all.fasta --out-file results/deeploc_partition_graphpart_needle_qc.csv 
python3 partitioning_quality/align_partitions.py --partition-file results/deeploc_partition_graphpart_mmseqs2.csv --fasta-file data/protein/DeepLoc/deeploc_all.fasta --out-file results/deeploc_partition_graphpart_mmseqs2_qc.csv


##
## 4. SignalP data
##
python3 baselines/homology_partition_cd_hit.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --threshold 0.3 --labels-name label --out-file results/signalp_partition_cd_hit.csv
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --threshold 0.3 --labels-name label --out-file results/signalp_reduction_cd_hit.csv

python3 baselines/homology_partition_mmseqs2.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --threshold 0.3 --labels-name label --out-file results/signalp_partition_mmseqs2_cmode0.csv --cluster-mode 0
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --threshold 0.3 --labels-name label --out-file results/signalp_reduction_mmseqs2_cmode0.csv --cluster-mode 0

graphpart needle -ff data/protein/SignalP5/train_set_reformatted.fasta -th 0.3 -ln label -of results/signalp_partition_graphpart_needle.csv -pa 5 --threads 20
graphpart mmseqs2 -ff data/protein/SignalP/train_set_reformatted.fasta -th 0.3 -ln label -of results/signalp_partition_graphpart_mmseqs2.csv -pa 5

python3 partitioning_quality/align_partitions.py --partition-file results/signalp_partition_cd_hit.csv --fasta-file data/protein/SignalP/train_set_reformatted.fasta --out-file results/signalp_partition_cd_hit_qc.csv --threads 20
python3 partitioning_quality/align_partitions.py --partition-file results/signalp_reduction_cd_hit.csv --fasta-file data/protein/SignalP/train_set_reformatted.fasta --out-file results/signalp_reduction_cd_hit_qc.csv --threads 20

python3 partitioning_quality/align_partitions.py --partition-file results/signalp_partition_mmseqs2_cmode0.csv --fasta-file data/protein/SignalP/train_set_reformatted.fasta --out-file results/signalp_partition_mmseqs2_cmode0_qc.csv --threads 15
python3 partitioning_quality/align_partitions.py --partition-file results/signalp_reduction_mmseqs2_cmode0.csv --fasta-file data/protein/SignalP/train_set_reformatted.fasta --out-file results/signalp_reduction_mmseqs2_cmode0_qc.csv --threads 15

python3 partitioning_quality/align_partitions.py --partition-file results/signalp_partition_graphpart_needle.csv --fasta-file data/protein/SignalP/train_set_reformatted.fasta --out-file results/signalp_partition_graphpart_needle_qc.csv 
python3 partitioning_quality/align_partitions.py --partition-file results/signalp_partition_graphpart_mmseqs2.csv --fasta-file data/protein/SignalP/train_set_reformatted.fasta --out-file results/signalp_partition_graphpart_mmseqs2_qc.csv


##
## 5. DeePromoter data
## 

# CD-HIT
python3 baselines/homology_partition_cd_hit.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --threshold 0.8 --out-file results/deepromoter_partition_cd_hit.csv --nucleotide
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --threshold 0.8 --out-file results/deepromoter_reduction_cd_hit.csv --nucleotide

python3 partitioning_quality/align_partitions.py --partition-file results/deepromoter_partition_cd_hit.csv --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --out-file results/deepromoter_partition_cd_hit_qc.csv --nucleotide --matrix EDNAFULL --threads 20
python3 partitioning_quality/align_partitions.py --partition-file results/deepromoter_reduction_cd_hit.csv --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --out-file results/deepromoter_reduction_cd_hit_qc.csv --nucleotide --matrix EDNAFULL --threads 20

# MMseqs2
python3 baselines/homology_partition_mmseqs2.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --threshold 0.8 --labels-name label --out-file results/deepromoter_partition_mmseqs2_cmode0.csv --cluster-mode 0
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --threshold 0.8 --labels-name label --out-file results/deepromoter_reduction_mmseqs2_cmode0.csv --cluster-mode 0

python3 partitioning_quality/align_partitions.py --partition-file results/deepromoter_partition_mmseqs2_cmode0.csv --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --out-file results/deepromoter_partition_mmseqs2_cmode0_qc.csv --nucleotide --matrix EDNAFULL --threads 15
python3 partitioning_quality/align_partitions.py --partition-file results/deepromoter_reduction_mmseqs2_cmode0.csv --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --out-file results/deepromoter_reduction_mmseqs2_cmode0_qc.csv --nucleotide --matrix EDNAFULL

# Graph-Part
graphpart needle -ff data/dna/DeePromoter/deepromoter_tata_human.fasta -th 0.8 -of results/deepromoter_partition_graphpart_needle.csv -pa 5 --threads 12 --nucleotide --matrix EDNAFULL
graphpart mmseqs2 -ff data/dna/DeePromoter/deepromoter_tata_human.fasta -th 0.8 -of results/deepromoter_partition_graphpart_mmseqs2.csv -pa 5 --nucleotide

python3 partitioning_quality/align_partitions.py --partition-file results/deepromoter_partition_graphpart_needle.csv --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --out-file results/deepromoter_partition_graphpart_needle_qc.csv --nucleotide --matrix EDNAFULL --threads 20
python3 partitioning_quality/align_partitions.py --partition-file results/deepromoter_partition_graphpart_mmseqs2.csv --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --out-file results/deepromoter_partition_graphpart_mmseqs2_qc.csv --nucleotide --matrix EDNAFULL


##
## 6. DeepM6ASeq data
## 

# CD-HIT
python3 baselines/homology_partition_cd_hit.py --fasta-file data/rna/DeepM6ASeq/deepm6aseq_dr.fasta --threshold 0.8 --out-file results/deepm6aseq_partition_cd_hit.csv --nucleotide --labels-name label
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/rna/DeepM6ASeq/deepm6aseq_dr.fasta --threshold 0.8 --out-file results/deepm6aseq_reduction_cd_hit.csv --nucleotide --labels-name label

python3 partitioning_quality/align_partitions.py --partition-file results/deepm6aseq_partition_cd_hit.csv --fasta-file data/rna/DeepM6ASeq/deepm6aseq_dr.fasta --out-file results/deepm6aseq_partition_cd_hit_qc.csv --nucleotide --matrix EDNAFULL --threads 20
python3 partitioning_quality/align_partitions.py --partition-file results/deepm6aseq_reduction_cd_hit.csv --fasta-file data/rna/DeepM6ASeq/deepm6aseq_dr.fasta --out-file results/deepm6aseq_reduction_cd_hit_qc.csv --nucleotide --matrix EDNAFULL --threads 20

# MMseqs2
python3 baselines/homology_partition_mmseqs2.py --fasta-file data/rna/DeepM6ASeq/deepm6aseq_dr.fasta --threshold 0.8 --labels-name label --out-file results/deepm6aseq_partition_mmseqs2_cmode0.csv --cluster-mode 0 --labels-name label
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/rna/DeepM6ASeq/deepm6aseq_dr.fasta --threshold 0.8 --labels-name label --out-file results/deepm6aseq_reduction_mmseqs2_cmode0.csv --cluster-mode 0 --labels-name label

python3 partitioning_quality/align_partitions.py --partition-file results/deepm6aseq_partition_mmseqs2_cmode0.csv --fasta-file data/rna/DeepM6ASeq/deepm6aseq_dr.fasta --out-file results/deepm6aseq_partition_mmseqs2_cmode0_qc.csv --nucleotide --matrix EDNAFULL --threads 15
python3 partitioning_quality/align_partitions.py --partition-file results/deepm6aseq_reduction_mmseqs2_cmode0.csv --fasta-file data/rna/DeepM6ASeq/deepm6aseq_dr.fasta --out-file results/deepm6aseq_reduction_mmseqs2_cmode0_qc.csv --nucleotide --matrix EDNAFULL

# Graph-Part
graphpart needle -ff data/rna/DeepM6ASeq/deepm6aseq_dr.fasta -th 0.8 -of results/deepm6aseq_partition_graphpart_needle.csv -pa 5 --threads 12 --nucleotide --matrix EDNAFULL --labels-name label
graphpart mmseqs2 -ff data/rna/DeepM6ASeq/deepm6aseq_dr.fasta -th 0.8 -of results/deepm6aseq_partition_graphpart_mmseqs2.csv -pa 5 --nucleotide --labels-name label

python3 partitioning_quality/align_partitions.py --partition-file results/deepm6aseq_partition_graphpart_needle.csv --fasta-file data/rna/DeepM6ASeq/deepm6aseq_dr.fasta --out-file results/deepm6aseq_partition_graphpart_needle_qc.csv --nucleotide --matrix EDNAFULL --threads 20
python3 partitioning_quality/align_partitions.py --partition-file results/deepm6aseq_partition_graphpart_mmseqs2.csv --fasta-file data/rna/DeepM6ASeq/deepm6aseq_dr.fasta --out-file results/deepm6aseq_partition_graphpart_mmseqs2_qc.csv --nucleotide --matrix EDNAFULL


##
## 7. Histone binding data
## 

# CD-HIT
python3 baselines/homology_partition_cd_hit.py --fasta-file data/dna/Histone/H3.fasta --threshold 0.8 --out-file results/histone_partition_cd_hit.csv --nucleotide --labels-name label
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/dna/Histone/H3.fasta --threshold 0.8 --out-file results/histone_reduction_cd_hit.csv --nucleotide --labels-name label

python3 partitioning_quality/align_partitions.py --partition-file results/histone_partition_cd_hit.csv --fasta-file data/dna/Histone/H3.fasta --out-file results/histone_partition_cd_hit_qc.csv --nucleotide --matrix EDNAFULL --threads 20
python3 partitioning_quality/align_partitions.py --partition-file results/histone_reduction_cd_hit.csv --fasta-file data/dna/Histone/H3.fasta --out-file results/histone_reduction_cd_hit_qc.csv --nucleotide --matrix EDNAFULL --threads 20

# MMseqs2
python3 baselines/homology_partition_mmseqs2.py --fasta-file data/dna/Histone/H3.fasta --threshold 0.8 --labels-name label --out-file results/histone_partition_mmseqs2_cmode0.csv --cluster-mode 0 --labels-name label
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/dna/Histone/H3.fasta --threshold 0.8 --labels-name label --out-file results/histone_reduction_mmseqs2_cmode0.csv --cluster-mode 0 --labels-name label

python3 partitioning_quality/align_partitions.py --partition-file results/histone_partition_mmseqs2_cmode0.csv --fasta-file data/dna/Histone/H3.fasta --out-file results/histone_partition_mmseqs2_cmode0_qc.csv --nucleotide --matrix EDNAFULL --threads 15
python3 partitioning_quality/align_partitions.py --partition-file results/histone_reduction_mmseqs2_cmode0.csv --fasta-file data/dna/Histone/H3.fasta --out-file results/histone_reduction_mmseqs2_cmode0_qc.csv --nucleotide --matrix EDNAFULL

# Graph-Part
graphpart needle -ff data/dna/Histone/H3.fasta -th 0.8 -of results/histone_partition_graphpart_needle.csv -pa 5 --threads 12 --nucleotide --matrix EDNAFULL --labels-name label
graphpart mmseqs2 -ff data/dna/Histone/H3.fasta -th 0.8 -of results/histone_partition_graphpart_mmseqs2.csv -pa 5 --nucleotide --labels-name label

python3 partitioning_quality/align_partitions.py --partition-file results/histone_partition_graphpart_needle.csv --fasta-file data/dna/Histone/H3.fasta --out-file results/histone_partition_graphpart_needle_qc.csv --nucleotide --matrix EDNAFULL --threads 20
python3 partitioning_quality/align_partitions.py --partition-file results/histone_partition_graphpart_mmseqs2.csv --fasta-file data/dna/Histone/H3.fasta --out-file results/histone_partition_graphpart_mmseqs2_qc.csv --nucleotide --matrix EDNAFULL


##
## 7. SPOT-RNA-1D data
## 

# CD-HIT
python3 baselines/homology_partition_cd_hit.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --threshold 0.8 --out-file results/spotrna1d_partition_cd_hit.csv --nucleotide --labels-name label
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --threshold 0.8 --out-file results/spotrna1d_reduction_cd_hit.csv --nucleotide --labels-name label

python3 partitioning_quality/align_partitions.py --partition-file results/spotrna1d_partition_cd_hit.csv --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --out-file results/spotrna1d_partition_cd_hit_qc.csv --nucleotide --matrix EDNAFULL --threads 20
python3 partitioning_quality/align_partitions.py --partition-file results/spotrna1d_reduction_cd_hit.csv --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --out-file results/spotrna1d_reduction_cd_hit_qc.csv --nucleotide --matrix EDNAFULL --threads 20

# MMseqs2
python3 baselines/homology_partition_mmseqs2.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --threshold 0.8 --labels-name label --out-file results/spotrna1d_partition_mmseqs2_cmode0.csv --cluster-mode 0 --labels-name label
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --threshold 0.8 --labels-name label --out-file results/spotrna1d_reduction_mmseqs2_cmode0.csv --cluster-mode 0 --labels-name label

python3 partitioning_quality/align_partitions.py --partition-file results/spotrna1d_partition_mmseqs2_cmode0.csv --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --out-file results/spotrna1d_partition_mmseqs2_cmode0_qc.csv --nucleotide --matrix EDNAFULL --threads 15
python3 partitioning_quality/align_partitions.py --partition-file results/spotrna1d_reduction_mmseqs2_cmode0.csv --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --out-file results/spotrna1d_reduction_mmseqs2_cmode0_qc.csv --nucleotide --matrix EDNAFULL

# Graph-Part
graphpart needle -ff data/rna/SPOT-RNA-1D/spot_rna_1d.fasta -th 0.8 -of results/spotrna1d_partition_graphpart_needle.csv -pa 5 --threads 12 --nucleotide --matrix EDNAFULL --labels-name label
graphpart mmseqs2 -ff data/rna/SPOT-RNA-1D/spot_rna_1d.fasta -th 0.8 -of results/spotrna1d_partition_graphpart_mmseqs2.csv -pa 5 --nucleotide --labels-name label

python3 partitioning_quality/align_partitions.py --partition-file results/spotrna1d_partition_graphpart_needle.csv --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --out-file results/spotrna1d_partition_graphpart_needle_qc.csv --nucleotide --matrix EDNAFULL --threads 20
python3 partitioning_quality/align_partitions.py --partition-file results/spotrna1d_partition_graphpart_mmseqs2.csv --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --out-file results/spotrna1d_partition_graphpart_mmseqs2_qc.csv --nucleotide --matrix EDNAFULL

