#!/bin/bash

cut -c1-50 ../data/protein/netgpi_dataset.fasta >50_aas.fasta
sed "s/>.*//" 50_aas.fasta > dummy_extension.fasta
paste -d "" 50_aas.fasta dummy_extension.fasta > 100_aas.fasta
paste -d "" 100_aas.fasta dummy_extension.fasta > 150_aas.fasta
paste -d "" 150_aas.fasta dummy_extension.fasta > 200_aas.fasta
paste -d "" 200_aas.fasta dummy_extension.fasta > 250_aas.fasta
paste -d "" 250_aas.fasta dummy_extension.fasta > 300_aas.fasta
paste -d "" 300_aas.fasta dummy_extension.fasta > 350_aas.fasta
paste -d "" 350_aas.fasta dummy_extension.fasta > 400_aas.fasta
paste -d "" 400_aas.fasta dummy_extension.fasta > 450_aas.fasta
paste -d "" 450_aas.fasta dummy_extension.fasta > 500_aas.fasta
paste -d "" 500_aas.fasta dummy_extension.fasta > 550_aas.fasta
paste -d "" 550_aas.fasta dummy_extension.fasta > 600_aas.fasta
paste -d "" 600_aas.fasta dummy_extension.fasta > 650_aas.fasta
paste -d "" 650_aas.fasta dummy_extension.fasta > 700_aas.fasta
paste -d "" 700_aas.fasta dummy_extension.fasta > 750_aas.fasta
paste -d "" 750_aas.fasta dummy_extension.fasta > 800_aas.fasta
paste -d "" 800_aas.fasta dummy_extension.fasta > 850_aas.fasta
paste -d "" 850_aas.fasta dummy_extension.fasta > 900_aas.fasta
paste -d "" 900_aas.fasta dummy_extension.fasta > 950_aas.fasta
paste -d "" 950_aas.fasta dummy_extension.fasta > 1000_aas.fasta



# Additional for mmseqs2
sed "s/>.*//" 1000_aas.fasta > dummy_extension.fasta

paste -d "" 1000_aas.fasta dummy_extension.fasta > 2000_aas.fasta
paste -d "" 2000_aas.fasta dummy_extension.fasta > 3000_aas.fasta
paste -d "" 3000_aas.fasta dummy_extension.fasta > 4000_aas.fasta
paste -d "" 4000_aas.fasta dummy_extension.fasta > 5000_aas.fasta
paste -d "" 5000_aas.fasta dummy_extension.fasta > 6000_aas.fasta
paste -d "" 6000_aas.fasta dummy_extension.fasta > 7000_aas.fasta
paste -d "" 7000_aas.fasta dummy_extension.fasta > 8000_aas.fasta
paste -d "" 8000_aas.fasta dummy_extension.fasta > 9000_aas.fasta
paste -d "" 9000_aas.fasta dummy_extension.fasta > 10000_aas.fasta