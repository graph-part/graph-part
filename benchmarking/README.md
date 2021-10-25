# Benchmarking experiments from manuscript

We benchmark the following aligmnent programs:  

- needleall
- CD-HIT / CD-HIT-EST
- MMseqs2

## Datasets

Dataset      | Type    | Number of classes | Orig. threshold | Download link 
-------------|---------|-------------------|-----------------|-----------------
NetGPI       | Protein | 2                 | 0.3             |
iLoc-mRNA    | RNA     | 5                 | 0.8             | TODO this is homology reduced already
DeepM6Aseq   | RNA     | 2                 | 0.8             | TODO this is homology reduced already
nRC          | RNA     | 13                | 0.8             | Script in `rna/`
DeepPromoter | DNA     | 2                 | None reported   | https://github.com/egochao/DeePromoter

### Data Processing

**nRC**  
Script handles everything.  
**DeepPromoter**  
Reformat to fasta  
`cat -n hs_pos_TATA.txt | sed 's/\ />seq_/' | sed 's/\t/|TATA\n/g' | sed 's/\ //g' >deepromoter_tata_human.fasta`


## TODO
DeepM6Aseq has sequence length of 102, too short?  
Alternative https://biodatamining.biomedcentral.com/articles/10.1186/s13040-017-0148-2#Sec2
also has more classes, so probably more interesting.  
Reimplement:

wget http://ftp.ebi.ac.uk/pub/databases/Rfam/14.6/fasta_files/RF00001.fa.gz




## Install additional dependencies

```
pip install scikit-learn
conda install -c conda-forge -c bioconda mmseqs2
conda install -c bioconda cd-hit
conda install -c bioconda blast
chmod + x baselines/psi_cd_hit.pl
chmod + x baselines/psi_cd_hit_local.pl
```


## Experiments

1. Graph-Part needleall partitioning
2. CD-HIT full cluster assignment
3. MMseqs2 full cluster assignment (different clustering modes)
4. MMseqs2 Graph-Part partitioning 



