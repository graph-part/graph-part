#!/bin/bash


# When upsampling, need to change identifier names.

for PREFIX in 'A' 'B' 'C' 'D' 'E' 'F'
do
sed "s/>/>${PREFIX}/" ../data/protein/netgpi_dataset.fasta >extension_${PREFIX}.fasta
done

cat extension* >source.fasta
cat source.fasta  | paste - - >source.tsv
rm extension*
rm source*

for NUMBER in 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000 13000 14000 15000 16000 17000 18000 19000 20000
do
shuf -n $NUMBER source.tsv | sed 's/\t/\n/' > "${NUMBER}_seqs.fasta"
done


