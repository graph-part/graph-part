# Benchmarking experiments from manuscript

We benchmark the following aligmnent programs:  

- needleall
- CD-HIT / CD-HIT-EST
- MMseqs2


This yields a total of six different experiments per dataset:

- Graph-Part needleall
- Graph-Part Mmseqs2
- Homology reduction CD-HIT
- Homology reduction MMseqs2
- Homology partitioning CD-HIT
- Homology reduction MMseqs2

## Datasets

Dataset      | Type    | Number of classes | Orig. threshold | Download link 
-------------|---------|-------------------|-----------------|-----------------
NetGPI       | Protein | 2                 | 0.3             | N/A
SignalP 5.0  | Protein | 4                 | 0.2             | https://services.healthtech.dtu.dk/service.php?SignalP-5.0  
DeepLoc 1.0  | Protein | 11                | 0.3             | https://services.healthtech.dtu.dk/service.php?DeepLoc-1.0
SPOT-RNA-1D  | RNA     | 1                 | 0.8             | Script in `data/rna/SPOT-RNA-1d`
DeepM6Aseq   | RNA     | 2                 | 0.8             | TODO this is homology reduced already
nRC          | RNA     | 13                | 0.8             | Script in `data/rna/nRC/`
DeepPromoter | DNA     | 1                 | None reported   | https://github.com/egochao/DeePromoter

### Data Processing

dataset directories have own readme files.


## Install additional dependencies

```
pip install scikit-learn
conda install -c conda-forge -c bioconda mmseqs2
conda install -c bioconda cd-hit
conda install -c bioconda blast
chmod + x baselines/psi_cd_hit.pl
chmod + x baselines/psi_cd_hit_local.pl
```


## Run benchmarking experiments

### Partitioning quality

We produce the results for partitioning quality in two steps. First, we compute partition assignments for all dataset using all algorithms, then we compute all cross-partition identities in order to find each sequence's maximum cross-partition identity (Figures 2,3). These results are also used to evaluate label balancing and runtime.

```bash
mkdir results
# this creates a .csv file for each dataset and algorithm
sh scripts/run_all_algorithms.sh
sh scripts/measure_cross_partition_identities.sh
```
### Runtime

The per-dataset runtimes are measured automatically using the `run_all_algorithms.sh` script. For the length and size dependency curves, we use the following commands:

```bash
mkdir runtime_benchmark
sh scripts/get_length_runtime_benchmarks_needle.sh
sh scripts/get_length_runtime_benchmarks_mmseqs2.sh

sh scripts/get_size_runtime_benchmarks_needle.sh
sh scripts/get_size_runtime_benchmarks_mmseqs2.sh
```

### Train-validation-test split comparison
This is the comparison of a 10-fold split with a 80%, 10%, 10% split, showing the utility of the train-validation-test split functionality to potentially retain more data.

```bash
mkdir results_splitmode
graphpart needle -ff data/protein/SignalP/train_set_reformatted.fasta -of results_splitmode/signalp_10fold.csv -ln label -th 0.3 -pa 10 -nt 24 -sc signalp.checkpoint
graphpart precomputed -ff data/protein/SignalP/train_set_reformatted.fasta -of results_splitmode/signalp_t01va01.csv -ln label -th 0.3 -ef signalp.checkpoint --test-ratio 0.1 --val-ratio 0.1
```
-----------------
### TODO delete everything below

### CD-HIT Reduction
```bash
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --out-file results/netgpi_reduction_cd_hit.csv --labels-name label --threshold 0.3
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --out-file results/deeploc_reduction_cd_hit.csv --labels-name label --threshold 0.3
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/protein/NetSurfP/netsurfp.fasta --out-file results/netsurfp_reduction_cd_hit.csv --threshold 0.25
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --out-file results/signalp_reduction_cd_hit.csv --labels-name label --threshold 0.3

python3 baselines/homology_reduce_cd_hit.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --out-file results/spotrna1d_reduction_cd_hit.csv --threshold 0.8 --nucleotide
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/rna/nRC/nRC_all.fasta --out-file results/nrc_reduction_cd_hit.csv --labels-name label --threshold 0.8 --nucleotide

python3 baselines/homology_reduce_cd_hit.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --out-file results/deepromoter_reduction_cd_hit.csv --threshold 0.8 --nucleotide
python3 baselines/homology_reduce_cd_hit.py --fasta-file data/dna/Histone/H3.fasta --out-file results/histone_reduction_cd_hit.csv --labels-name label --threshold 0.8 --nucleotide
```

### CD-HIT Partitioning
```bash
python3 baselines/homology_partition_cd_hit.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --out-file results/netgpi_partition_cd_hit.csv --labels-name label --threshold 0.3
python3 baselines/homology_partition_cd_hit.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --out-file results/deeploc_partition_cd_hit.csv --labels-name label --threshold 0.3
python3 baselines/homology_partition_cd_hit.py --fasta-file data/protein/NetSurfP/netsurfp.fasta --out-file results/netsurfp_partition_cd_hit.csv --threshold 0.25
python3 baselines/homology_partition_cd_hit.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --out-file results/signalp_partition_cd_hit.csv --labels-name label --threshold 0.3

python3 baselines/homology_partition_cd_hit.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --out-file results/spotrna1d_partition_cd_hit.csv --threshold 0.8 --nucleotide
python3 baselines/homology_partition_cd_hit.py --fasta-file data/rna/nRC/nRC_all.fasta --out-file results/nrc_partition_cd_hit.csv --labels-name label --threshold 0.8 --nucleotide

python3 baselines/homology_partition_cd_hit.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --out-file results/deepromoter_partition_cd_hit.csv --threshold 0.8 --nucleotide
python3 baselines/homology_partition_cd_hit.py --fasta-file data/dna/Histone/H3.fasta --out-file results/histone_partition_cd_hit.csv --labels-name label --threshold 0.8 --nucleotide
```

### MMseqs2 Reduction
```bash
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --out-file results/netgpi_reduction_mmseqs2_cmode0.csv --labels-name label --threshold 0.3
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --out-file results/deeploc_reduction_mmseqs2_cmode0.csv --labels-name label --threshold 0.3
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/protein/NetSurfP/netsurfp.fasta --out-file results/netsurfp_reduction_mmseqs2_cmode0.csv --threshold 0.25
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --out-file results/signalp_reduction_mmseqs2_cmode0.csv --labels-name label --threshold 0.3

python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --out-file results/spotrna1d_reduction_mmseqs2_cmode0.csv --threshold 0.8 --nucleotide
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/rna/nRC/nRC_all.fasta --out-file results/nrc_reduction_mmseqs2_cmode0.csv --labels-name label --threshold 0.8 --nucleotide

python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --out-file results/deepromoter_reduction_mmseqs2_cmode0.csv --threshold 0.8 --nucleotide
python3 baselines/homology_reduce_mmseqs2.py --fasta-file data/dna/Histone/H3.fasta --out-file results/histone_reduction_mmseqs2_cmode0.csv --labels-name label --threshold 0.8 --nucleotide
```

### MMseqs2 Partitioning
```bash
python3 baselines/homology_partition_mmseqs2.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --out-file results/netgpi_partition_mmseqs2_cmode0.csv --labels-name label --threshold 0.3
python3 baselines/homology_partition_mmseqs2.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --out-file results/deeploc_partition_mmseqs2_cmode0.csv --labels-name label --threshold 0.3
python3 baselines/homology_partition_mmseqs2.py --fasta-file data/protein/NetSurfP/netsurfp.fasta --out-file results/netsurfp_partition_mmseqs2_cmode0.csv --threshold 0.25
python3 baselines/homology_partition_mmseqs2.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --out-file results/signalp_partition_mmseqs2_cmode0.csv --labels-name label --threshold 0.3

python3 baselines/homology_partition_mmseqs2.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --out-file results/spotrna1d_partition_mmseqs2_cmode0.csv --threshold 0.8 --nucleotide
python3 baselines/homology_partition_mmseqs2.py --fasta-file data/rna/nRC/nRC_all.fasta --out-file results/nrc_partition_mmseqs2_cmode0.csv --labels-name label --threshold 0.8 --nucleotide

python3 baselines/homology_partition_mmseqs2.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --out-file results/deepromoter_partition_mmseqs2_cmode0.csv --threshold 0.8 --nucleotide
python3 baselines/homology_partition_mmseqs2.py --fasta-file data/dna/Histone/H3.fasta --out-file results/histone_partition_mmseqs2_cmode0.csv --labels-name label --threshold 0.8 --nucleotide
```

### GraphPart needle
```bash
graphpart needle -ff data/protein/NetGPI/netgpi_dataset.fasta -of results/netgpi_partition_graphpart_needle.csv -ln label -th 0.3 -nt 64
graphpart needle -ff data/protein/DeepLoc/deeploc_all.fasta -of results/deeploc_partition_graphpart_needle.csv -ln label -th 0.3 -nt 64 -sc deeploc.checkpoint
graphpart needle -ff data/protein/NetSurfP/netsurfp.fasta -of results/netsurfp_partition_graphpart_needle.csv -th 0.25 -nt 64
graphpart needle -ff data/protein/SignalP/train_set_reformatted.fasta -of results/signalp_partition_graphpart_needle.csv -ln label -th 0.3 -nt 64

graphpart needle -ff data/rna/SPOT-RNA-1D/spot_rna_1d.fasta -of results/spotrna1d_partition_graphpart_needle.csv -th 0.8 -nu -nt 64
graphpart needle -ff data/rna/nRC/nRC_all.fasta -of results/nrc_partition_graphpart_needle.csv -ln label -th 0.8 -nu -nt 64

graphpart needle -ff data/dna/DeePromoter/deepromoter_tata_human.fasta -of results/deepromoter_partition_graphpart_needle.csv -th 0.8 -nu -nt 64
graphpart needle -ff data/dna/Histone/H3.fasta -of results/histone_partition_graphpart_needle.csv -ln label -th 0.8 -nu -nt 64
```

### GraphPart mmseqs2
```bash
graphpart mmseqs2 -ff data/protein/NetGPI/netgpi_dataset.fasta -of results/netgpi_partition_graphpart_mmseqs2.csv -ln label -th 0.3 -nt 64
graphpart mmseqs2 -ff data/protein/DeepLoc/deeploc_all.fasta -of results/deeploc_partition_graphpart_mmseqs2.csv -ln label -th 0.3 -nt 64 -sc deeploc.checkpoint
graphpart mmseqs2 -ff data/protein/NetSurfP/netsurfp.fasta -of results/netsurfp_partition_graphpart_mmseqs2.csv -th 0.25 -nt 64
graphpart mmseqs2 -ff data/protein/SignalP/train_set_reformatted.fasta -of results/signalp_partition_graphpart_mmseqs2.csv -ln label -th 0.3 -nt 64

graphpart mmseqs2 -ff data/rna/SPOT-RNA-1D/spot_rna_1d.fasta -of results/spotrna1d_partition_graphpart_mmseqs2.csv -th 0.8 -nu -nt 64
graphpart mmseqs2 -ff data/rna/nRC/nRC_all.fasta -of results/nrc_partition_graphpart_mmseqs2.csv -ln label -th 0.8 -nu -nt 64

graphpart mmseqs2 -ff data/dna/DeePromoter/deepromoter_tata_human.fasta -of results/deepromoter_partition_graphpart_mmseqs2.csv -th 0.8 -nu -nt 64
graphpart mmseqs2 -ff data/dna/Histone/H3.fasta -of results/histone_partition_graphpart_mmseqs2.csv -ln label -th 0.8 -nu -nt 64
```
## GraphPart mmseqs2needle
```bash
graphpart mmseqs2needle -ff data/protein/NetGPI/netgpi_dataset.fasta -of results_mmseqs2needle/netgpi_t_01.csv -ln label -th 0.3 -re 0.4 -nt 64 -re 0.1 # Recomputing exact values for 683520 pairs.
```

## sklearn
```bash
python3 baselines/random_partition_sklearn.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --out-file results/netgpi_sklearn.csv --labels-name label
python3 baselines/random_partition_sklearn.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --out-file results/deeploc_sklearn.csv --labels-name label
python3 baselines/random_partition_sklearn.py --fasta-file data/protein/NetSurfP/netsurfp.fasta --out-file results/netsurfp_sklearn.csv
python3 baselines/random_partition_sklearn.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --out-file results/signalp_sklearn.csv --labels-name label

python3 baselines/random_partition_sklearn.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --out-file results/spotrna1d_sklearn.csv --nucleotide
python3 baselines/random_partition_sklearn.py --fasta-file data/rna/nRC/nRC_all.fasta --out-file results/nrc_sklearn.csv --labels-name label --nucleotide

python3 baselines/random_partition_sklearn.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --out-file results/deepromoter_sklearn.csv --nucleotide
python3 baselines/random_partition_sklearn.py --fasta-file data/dna/Histone/H3.fasta --out-file results/histone_sklearn.csv --labels-name label --nucleotide
```

### Evaluate
```bash

python3 partitioning_quality/align_partitions.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --partition-file results/netgpi_partition_graphpart_needle.csv --out-file results/netgpi_partition_graphpart_neeedle_qc.csv


python3 partitioning_quality/align_partitions.py --fasta-file data/protein/NetGPI/netgpi_dataset.fasta --partition-file results/netgpi_sklearn.csv --out-file results/netgpi_sklearn_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/DeepLoc/deeploc_all.fasta --partition-file results/deeploc_sklearn.csv --out-file results/deeploc_sklearn_qc.csv --threads 20
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/NetSurfP/netsurfp.fasta --partition-file results/netsurfp_sklearn.csv --out-file results/netsurfp_sklearn_qc.csv
python3 partitioning_quality/align_partitions.py --fasta-file data/protein/SignalP/train_set_reformatted.fasta --partition-file results/signalp_sklearn.csv --out-file results/signalp_sklearn_qc.csv

python3 partitioning_quality/align_partitions.py --fasta-file data/rna/SPOT-RNA-1D/spot_rna_1d.fasta --partition-file results/spotrna1d_sklearn.csv --out-file results/spotrna1d_sklearn_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/rna/nRC/nRC_all.fasta --partition-file results/nrc_sklearn.csv --out-file results/nrc_sklearn_qc.csv --nucleotide --matrix EDNAFULL

python3 partitioning_quality/align_partitions.py --fasta-file data/dna/DeePromoter/deepromoter_tata_human.fasta --partition-file results/deepromoter_sklearn.csv --out-file results/deepromoter_sklearn_qc.csv --nucleotide --matrix EDNAFULL
python3 partitioning_quality/align_partitions.py --fasta-file data/dna/Histone/H3.fasta --partition-file results/histone_sklearn.csv --out-file results/histone_sklearn_qc.csv --nucleotide --matrix EDNAFULL
```


### Train-validation-test split comparison
```bash
graphpart precomputed -ff data/protein/DeepLoc/deeploc_all.fasta -of results_splitmode/deeploc_10fold.csv -ln label -th 0.3 -ef deeploc.checkpoint -pa 10
graphpart precomputed -ff data/protein/DeepLoc/deeploc_all.fasta -of results_splitmode/deeploc_t01va01.csv -ln label -th 0.3 -ef deeploc.checkpoint --test-ratio 0.1 --val-ratio 0.1

graphpart needle -ff data/protein/NetGPI/netgpi_dataset.fasta -of results_splitmode/netgpi_10fold.csv -ln label -th 0.3 -pa 10 -nt 24 -sc netgpi.checkpoint
graphpart precomputed -ff data/protein/NetGPI/netgpi_dataset.fasta -of results_splitmode/netgpi_10fold_nm.csv -ln label -th 0.3 -ef netgpi.checkpoint -pa 10 -nm
graphpart precomputed -ff data/protein/NetGPI/netgpi_dataset.fasta -of results_splitmode/netgpi_t01va01.csv -ln label -th 0.3 -ef netgpi.checkpoint --test-ratio 0.1 --val-ratio 0.1

graphpart needle -ff data/protein/SignalP/train_set_reformatted.fasta -of results_splitmode/signalp_10fold.csv -ln label -th 0.3 -pa 10 -nt 24 -sc signalp.checkpoint
graphpart precomputed -ff data/protein/SignalP/train_set_reformatted.fasta -of results_splitmode/signalp_10fold_nm.csv -ln label -th 0.3 -ef signalp.checkpoint -pa 10 -nm
graphpart precomputed -ff data/protein/SignalP/train_set_reformatted.fasta -of results_splitmode/signalp_t01va01.csv -ln label -th 0.3 -ef signalp.checkpoint --test-ratio 0.1 --val-ratio 0.1
graphpart precomputed -ff data/protein/SignalP/train_set_reformatted.fasta -of results_splitmode/signalp_t01va01_nm.csv -ln label -th 0.3 -ef signalp.checkpoint --test-ratio 0.1 --val-ratio 0.1 -nm

graphpart needle -ff data/protein/NetSurfP/netsurfp.fasta -of results_splitmode/netsurfp_10fold.csv -th 0.25 -pa 10 -nt 24 -sc netsurfp.checkpoint


```

### find t'
graphpart mmseqs2needle -ff data/protein/NetGPI/netgpi_dataset.fasta -of results_mmseqs2needle/netgpi_t_01.csv -ln label -th 0.3 -re 0.4 -nt 32 -re 0.1 -pm multiprocess