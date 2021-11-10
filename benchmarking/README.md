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






