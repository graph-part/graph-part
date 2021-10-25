#!/bin/bash


# When upsampling, need to change identifier names.
# in each loop, prefix previous identifier with another letter.

for PREFIX in 'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L'
do
sed "s/>/>${PREFIX}/" ../data/protein/netgpi_dataset.fasta >extension_${PREFIX}.fasta
done

cat extension* >source.fasta
cat source.fasta  | paste - - >source.tsv
rm extension*

for NUMBER in 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000 13000 14000 15000 16000 17000 18000 19000 20000 21000 22000 23000 24000 25000 26000 27000 28000 29000 30000
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
cat source.fasta  | paste - - >source.tsv


for NUMBER in 35000 40000 45000 50000 55000 65000 70000 75000 80000 85000 90000 95000 100000
do
shuf -n $NUMBER source.tsv | sed 's/\t/\n/' > "${NUMBER}_seqs.fasta"
done

rm extension*
rm source*