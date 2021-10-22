#!/bin/bash


# When upsampling, need to change identifier names.
# in each loop, prefix previous identifier with another letter.

for PREFIX in 'A' 'B' 'C' 'D' 'E' 'F'
do
sed "s/>/>${PREFIX}/" ../data/protein/netgpi_dataset.fasta >extension_${PREFIX}.fasta
done

cat extension* >source.fasta
cat source.fasta  | paste - - >source.tsv
rm extension*

for NUMBER in 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000 13000 14000 15000 16000 17000 18000 19000 20000
do
shuf -n $NUMBER source.tsv | sed 's/\t/\n/' > "${NUMBER}_seqs.fasta"
done



# Additional for mmseqs2

# this goes to 120k
for PREFIX in 'A' 'B' 'C' 'D' 'E' 'F'
do
sed "s/>/>${PREFIX}/" source.fasta >extension_${PREFIX}.fasta
done

cat extension* >source.fasta

# this brings us over 500k
for PREFIX in 'A' 'B' 'C' 'D' 'E' 
do
sed "s/>/>${PREFIX}/" source.fasta >extension_${PREFIX}.fasta
done

cat extension* >source.fasta


cat source.fasta  | paste - - >source.tsv
rm extension*

for NUMBER in 50000 100000 150000 200000 250000 300000 450000 500000
do
shuf -n $NUMBER source.tsv | sed 's/\t/\n/' > "${NUMBER}_seqs.fasta"
done


rm source*