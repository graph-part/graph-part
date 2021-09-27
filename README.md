# Graph-Part
Protein dataset partitioning pipeline (Gíslason 2021)


# Instructions

## Input format
The fasta headers contain the metadata used for clustering. Graph-Part expects the following encoding:
```
>P35583|label=EUKARYANO_SP
```
We're not using the priority labeling function at the moment.


### Step 1
Generate pairwise identities using ggsearch36
```
python3 get_edgelist.py --sequences input.fasta --ggs '/work3/felteu/fasta36/bin/ggsearch36' --outfile edgelist.csv
```

### Step 2
Use metadata from fasta headers + pairwise identities to obtain clustering.

```
python3 graphpart.py --meta-file input.fasta --edge-file edgelist.csv --threshold 0.3 --transformation one-minus --out-file graphpart_assignments.csv \
 --labels-name label --partitions 3
```