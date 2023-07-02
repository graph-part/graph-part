# Benchmarking experiments from manuscript

We benchmark the following methods:  

- GraphPart needle / mmseqs2
- PSI-CD-HIT / CD-HIT-EST
- MMseqs2
- Hobohm 2 needle / mmseqs2


This yields a total of eight different experiments per dataset:

- Graph-Part needleall
- Graph-Part Mmseqs2
- Homology reduction CD-HIT
- Homology reduction MMseqs2
- Homology partitioning CD-HIT
- Homology reduction MMseqs2
- Homology reduction Hobohm 2 needleall
- Homology reduction Hobohm 2 MMseqs2

## Datasets

Dataset      | Type    | Number of classes | Orig. threshold | Download link 
-------------|---------|-------------------|-----------------|-----------------
NetGPI       | Protein | 2                 | 0.3             | https://services.healthtech.dtu.dk/services/NetGPI-1.1/
SignalP 5.0  | Protein | 4                 | 0.2             | https://services.healthtech.dtu.dk/service.php?SignalP-5.0  
DeepLoc 1.0  | Protein | 11                | 0.3             | https://services.healthtech.dtu.dk/service.php?DeepLoc-1.0
NetSurfP 2.1 | Protein | 1                 | 0.25            | Described in [`data/protein/NetSurfP`](data/protein/NetSurfP/README.md)
SPOT-RNA-1D  | RNA     | 1                 | 0.8             | Described in [`data/rna/SPOT-RNA-1d`](data/rna/SPOT-RNA-1D/README.md)
nRC          | RNA     | 13                | 0.8             | Described in [`data/rna/nRC/`](data/rna/nRC/README.md)
DeepPromoter | DNA     | 1                 | None reported   | https://github.com/egochao/DeePromoter
Histone H3   | DNA     | 2                 | None reported   | https://github.com/Doulrs/Hilbert-CNN

### Data Processing

The dataset directories have README files.


## Install additional dependencies

```sh
pip install scikit-learn
conda install -c conda-forge -c bioconda mmseqs2
conda install -c bioconda cd-hit
conda install -c bioconda blast
chmod +x baselines/psi_cd_hit.pl
chmod +x baselines/psi_cd_hit_local.pl
```


## Run benchmarking experiments

### Partitioning quality

We produce the results for partitioning quality in two steps. First, we compute partition assignments for all datasets using all algorithms, then we compute all cross-partition identities in order to find each sequence's maximum cross-partition identity (Figures 2,3). These results are also used to evaluate label balancing and runtime.

```bash
mkdir results
# this creates a .csv file for each dataset and algorithm
sh scripts/run_all_algorithms.sh
# this creates a _qc.csv file for each dataset and algorithm.
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
# with sequence moving
graphpart precomputed -ff data/protein/SignalP/train_set_reformatted.fasta -of results_splitmode/signalp_t01va01.csv -ln label -th 0.3 -ef signalp.checkpoint --test-ratio 0.1 --val-ratio 0.1
# without
graphpart precomputed -ff data/protein/SignalP/train_set_reformatted.fasta -of results_splitmode/signalp_t01va01_nm.csv -ln label -th 0.3 -ef signalp.checkpoint --test-ratio 0.1 --val-ratio 0.1 -nm
```