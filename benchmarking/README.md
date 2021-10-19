# Benchmarking experiments from manuscript

We benchmark the following aligmnent programs:  

- needleall
- CD-HIT / CD-HIT-EST
- MMseqs2

## Datasets

Dataset   | Type    | Number of classes | Threshold | Link 
----------|---------|-------------------|-----------|-----------------
NetGPI    | Protein | 2                 | 0.3       |
iLoc-mRNA | RNA     | 5                 | 0.8       |
 |DNA|

## Install additional dependencies
```
conda install -c conda-forge -c bioconda mmseqs2
conda install -c bioconda cd-hit
```


## Experiments

1. Ass