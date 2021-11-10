PDB Query:
```

QUERY: Number of Protein Instances (Chains) per Assembly < 1 AND Polymer Entity Type = "RNA" AND Refinement Resolution <= 3.5 AND Entry Polymer Composition = "RNA"
```

Download from PDB:  

```
./batch_download.sh -f list_file.txt -p
mkdir raw
mv *.gz raw/
gunzip raw/*
```



Deprecated, just use Python script instead.
```
pip install rna-tools
# split chains and get one fasta for all
rna_pdb_toolsx.py --get-seq raw/* --fasta >pdb_all.fasta

# make tab separated
# fix upper lower case in seqs
# drop protein sequences
awk '{if($2~/[AUGTCaugtc]*/) {print toupper($1)"\t"$2}}'
# drop heteroatoms at start and end
# drop sequences with gaps or heteroatoms in between

# awk one-line magic
#drops non-nucleotide seqs 
#fixes to upper case.
#removes path from the fasta header
awk -v RS=">" -v FS="\n" '{if($2~/^[AUGTCaugtc]*$/) {n=split($1, header, "/"); print ">"header[n]"\n"toupper($2)}}'
```
