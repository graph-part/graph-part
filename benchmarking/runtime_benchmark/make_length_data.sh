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