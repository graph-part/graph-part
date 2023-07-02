
# GraphPart
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --partition-file results/netgpi_partition_graphpart_needle.csv --out-file results/netgpi_partition_graphpart_neeedle_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --partition-file results/netgpi_partition_graphpart_mmseqs2.csv --out-file results/netgpi_partition_graphpart_mmseqs2_qc.csv

python3 partitioning_quality/align_partitions.py --fasta-file data/rna/nRC/nRC_all.fasta --partition-file results/nrc_partition_graphpart_needle.csv --out-file results/nrc_partition_graphpart_neeedle_qc.csv  --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/rna/nRC/nRC_all.fasta --partition-file results/nrc_partition_graphpart_mmseqs2.csv --out-file results/nrc_partition_graphpart_mmseqs2_qc.csv --nucleotide --matrix EDNAFULL

python3 partitioning_quality/align_partitions.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --partition-file results/deeploc_partition_graphpart_needle.csv --out-file results/deeploc_partition_graphpart_neeedle_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --partition-file results/deeploc_partition_graphpart_mmseqs2.csv --out-file results/deeploc_partition_graphpart_mmseqs2_qc.csv

python3 partitioning_quality/align_partitions.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --partition-file results/signalp_partition_graphpart_needle.csv --out-file results/signalp_partition_graphpart_neeedle_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --partition-file results/signalp_partition_graphpart_mmseqs2.csv --out-file results/signalp_partition_graphpart_mmseqs2_qc.csv

python3 partitioning_quality/align_partitions.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --partition-file results/deepromoter_partition_graphpart_needle.csv --out-file results/deepromoter_partition_graphpart_neeedle_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --partition-file results/deepromoter_partition_graphpart_mmseqs2.csv --out-file results/deepromoter_partition_graphpart_mmseqs2_qc.csv --nucleotide --matrix EDNAFULL

python3 partitioning_quality/align_partitions.py --fasta-file data/dna/Histone/H3.fasta --partition-file results/histone_partition_graphpart_needle.csv --out-file results/histone_partition_graphpart_neeedle_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/dna/Histone/H3.fasta --partition-file results/histone_partition_graphpart_mmseqs2.csv --out-file results/histone_partition_graphpart_mmseqs2_qc.csv --nucleotide --matrix EDNAFULL

python3 partitioning_quality/align_partitions.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --partition-file results/spotrna1d_partition_graphpart_needle.csv --out-file results/spotrna1d_partition_graphpart_neeedle_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --partition-file results/spotrna1d_partition_graphpart_mmseqs2.csv --out-file results/spotrna1d_partition_graphpart_mmseqs2_qc.csv --nucleotide --matrix EDNAFULL

python3 partitioning_quality/align_partitions.py --fasta-file data/protein/NetSurfP/netsurfp.fasta --partition-file results/netsurfp_partition_graphpart_needle.csv --out-file results/netsurfp_partition_graphpart_neeedle_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/NetSurfP/netsurfp.fasta --partition-file results/netsurfp_partition_graphpart_mmseqs2.csv --out-file results/netsurfp_partition_graphpart_mmseqs2_qc.csv

# CD-HIT partitioning
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --partition-file results/netgpi_partition_cdhit.csv --out-file results/netgpi_partition_cdhit_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/rna/nRC/nRC_all.fasta --partition-file results/nrc_partition_cdhit.csv --out-file results/nrc_partition_cdhit_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --partition-file results/deeploc_partition_cdhit.csv --out-file results/deeploc_partition_cdhit_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --partition-file results/signalp_partition_cdhit.csv --out-file results/signalp_partition_cdhit_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --partition-file results/deepromoter_partition_cdhit.csv --out-file results/deepromoter_partition_cdhit_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/dna/Histone/H3.fasta --partition-file results/histone_partition_cdhit.csv --out-file results/histone_partition_cdhit_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --partition-file results/spotrna1d_partition_cdhit.csv --out-file results/spotrna1d_partition_cdhit_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/NetSurfP/netsurfp.fasta --partition-file results/netsurfp_partition_cdhit.csv --out-file results/netsurfp_partition_cdhit_qc.csv

# MMseqs2 partitioning
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --partition-file results/netgpi_partition_mmseqs2.csv --out-file results/netgpi_partition_mmseqs2_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/rna/nRC/nRC_all.fasta --partition-file results/nrc_partition_mmseqs2.csv --out-file results/nrc_partition_mmseqs2_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --partition-file results/deeploc_partition_mmseqs2.csv --out-file results/deeploc_partition_mmseqs2_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --partition-file results/signalp_partition_mmseqs2.csv --out-file results/signalp_partition_mmseqs2_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --partition-file results/deepromoter_partition_mmseqs2.csv --out-file results/deepromoter_partition_mmseqs2_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/dna/Histone/H3.fasta --partition-file results/histone_partition_mmseqs2.csv --out-file results/histone_partition_mmseqs2_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --partition-file results/spotrna1d_partition_mmseqs2.csv --out-file results/spotrna1d_partition_mmseqs2_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/NetSurfP/netsurfp.fasta --partition-file results/netsurfp_partition_mmseqs2.csv --out-file results/netsurfp_partition_mmseqs2_qc.csv

# CD-HIT reduction
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --partition-file results/netgpi_reduction_cdhit.csv --out-file results/netgpi_reduction_cdhit_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/rna/nRC/nRC_all.fasta --partition-file results/nrc_reduction_cdhit.csv --out-file results/nrc_reduction_cdhit_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --partition-file results/deeploc_reduction_cdhit.csv --out-file results/deeploc_reduction_cdhit_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --partition-file results/signalp_reduction_cdhit.csv --out-file results/signalp_reduction_cdhit_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --partition-file results/deepromoter_reduction_cdhit.csv --out-file results/deepromoter_reduction_cdhit_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/dna/Histone/H3.fasta --partition-file results/histone_reduction_cdhit.csv --out-file results/histone_reduction_cdhit_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --partition-file results/spotrna1d_reduction_cdhit.csv --out-file results/spotrna1d_reduction_cdhit_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/NetSurfP/netsurfp.fasta --partition-file results/netsurfp_reduction_cdhit.csv --out-file results/netsurfp_reduction_cdhit_qc.csv

# MMseqs2 reduction
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --partition-file results/netgpi_reduction_mmseqs2.csv --out-file results/netgpi_reduction_mmseqs2_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/rna/nRC/nRC_all.fasta --partition-file results/nrc_reduction_mmseqs2.csv --out-file results/nrc_reduction_mmseqs2_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --partition-file results/deeploc_reduction_mmseqs2.csv --out-file results/deeploc_reduction_mmseqs2_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --partition-file results/signalp_reduction_mmseqs2.csv --out-file results/signalp_reduction_mmseqs2_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --partition-file results/deepromoter_reduction_mmseqs2.csv --out-file results/deepromoter_reduction_mmseqs2_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/dna/Histone/H3.fasta --partition-file results/histone_reduction_mmseqs2.csv --out-file results/histone_reduction_mmseqs2_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --partition-file results/spotrna1d_reduction_mmseqs2.csv --out-file results/spotrna1d_reduction_mmseqs2_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/NetSurfP/netsurfp.fasta --partition-file results/netsurfp_reduction_mmseqs2.csv --out-file results/netsurfp_reduction_mmseqs2_qc.csv

# sklearn partitioning
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --partition-file results/netgpi_sklearn.csv --out-file results/netgpi_sklearn_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/rna/nRC/nRC_all.fasta --partition-file results/nrc_sklearn.csv --out-file results/nrc_sklearn_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --partition-file results/deeploc_sklearn.csv --out-file results/deeploc_sklearn_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --partition-file results/signalp_sklearn.csv --out-file results/signalp_sklearn_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --partition-file results/deepromoter_sklearn.csv --out-file results/deepromoter_sklearn_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/dna/Histone/H3.fasta --partition-file results/histone_sklearn.csv --out-file results/histone_sklearn_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --partition-file results/spotrna1d_sklearn.csv --out-file results/spotrna1d_sklearn_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/NetSurfP/netsurfp.fasta --partition-file results/netsurfp_sklearn.csv --out-file results/netsurfp_sklearn_qc.csv

python3 partitioning_quality/align_partitions.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --partition-file results/netgpi_reduction_hobohm2.csv --out-file results/netgpi_reduction_hobohm2_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/rna/nRC/nRC_all.fasta --partition-file results/nrc_reduction_hobohm2.csv --out-file results/nrc_reduction_hobohm2_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --partition-file results/deeploc_reduction_hobohm2.csv --out-file results/deeploc_reduction_hobohm2_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --partition-file results/signalp_reduction_hobohm2.csv --out-file results/signalp_reduction_hobohm2_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --partition-file results/deepromoter_reduction_hobohm2.csv --out-file results/deepromoter_reduction_hobohm2_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/dna/Histone/H3.fasta --partition-file results/histone_reduction_hobohm2.csv --out-file results/histone_reduction_hobohm2_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --partition-file results/spotrna1d_reduction_hobohm2.csv --out-file results/spotrna1d_reduction_hobohm2_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/NetSurfP/netsurfp.fasta --partition-file results/netsurfp_reduction_hobohm2.csv --out-file results/netsurfp_reduction_hobohm2_qc.csv


# Hobohm 2 reduction

python3 partitioning_quality/align_partitions.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --partition-file results/netgpi_reduction_hobohm2.csv --out-file results/netgpi_reduction_hobohm2_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/rna/nRC/nRC_all.fasta --partition-file results/nrc_reduction_hobohm2.csv --out-file results/nrc_reduction_hobohm2_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --partition-file results/deeploc_reduction_hobohm2.csv --out-file results/deeploc_reduction_hobohm2_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --partition-file results/signalp_reduction_hobohm2.csv --out-file results/signalp_reduction_hobohm2_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --partition-file results/deepromoter_reduction_hobohm2.csv --out-file results/deepromoter_reduction_hobohm2_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/dna/Histone/H3.fasta --partition-file results/histone_reduction_hobohm2.csv --out-file results/histone_reduction_hobohm2_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --partition-file results/spotrna1d_reduction_hobohm2.csv --out-file results/spotrna1d_reduction_hobohm2_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/NetSurfP/netsurfp.fasta --partition-file results/netsurfp_reduction_hobohm2.csv --out-file results/netsurfp_reduction_hobohm2_qc.csv

python3 partitioning_quality/align_partitions.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --partition-file results/netgpi_reduction_hobohm2_mmseqs2.csv --out-file results/netgpi_reduction_hobohm2_mmseqs2_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/rna/nRC/nRC_all.fasta --partition-file results/nrc_reduction_hobohm2_mmseqs2.csv --out-file results/nrc_reduction_hobohm2_mmseqs2_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --partition-file results/deeploc_reduction_hobohm2_mmseqs2.csv --out-file results/deeploc_reduction_hobohm2_mmseqs2_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --partition-file results/signalp_reduction_hobohm2_mmseqs2.csv --out-file results/signalp_reduction_hobohm2_mmseqs2_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --partition-file results/deepromoter_reduction_hobohm2_mmseqs2.csv --out-file results/deepromoter_reduction_hobohm2_mmseqs2_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/dna/Histone/H3.fasta --partition-file results/histone_reduction_hobohm2_mmseqs2.csv --out-file results/histone_reduction_hobohm2_mmseqs2_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --partition-file results/spotrna1d_reduction_hobohm2_mmseqs2.csv --out-file results/spotrna1d_reduction_hobohm2_mmseqs2_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/NetSurfP/netsurfp.fasta --partition-file results/netsurfp_reduction_hobohm2_mmseqs2.csv --out-file results/netsurfp_reduction_hobohm2_mmseqs2_qc.csv
