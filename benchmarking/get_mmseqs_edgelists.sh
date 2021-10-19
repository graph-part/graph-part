#!/bin/bash

conda activate gp-env

mkdir temp
# 1. NetGPI
mmseqs createdb data/netgpi_dataset.fasta temp/netgpi_db
mmseqs prefilter -s 7.5 temp/netgpi_db temp/netgpi_db temp/netgpi_pref
mmseqs align temp/netgpi_db temp/netgpi_db temp/netgpi_pref temp/netgpi_align_db
mmseqs convertalis temp/netgpi_db temp/netgpi_db temp/netgpi_align_db temp/alignments.tab # column 3 contains sequence identity.
awk -F'[|\t]' '{print $1","$4","$7}' temp/alignments.tab >netgpi_mmseqs_edgelist.csv


# 2. iLoc-mRNA

mmseqs createdb --dbtype 2 data/mloci_dataset.fasta temp/mloci_db
mmseqs prefilter -s 7.5 temp/mloci_db temp/mloci_db temp/mloci_pref
mmseqs align temp/mloci_db temp/mloci_db temp/mloci_pref temp/mloci_align_db
mmseqs convertalis temp/mloci_db temp/mloci_db temp/mloci_align_db temp/alignments.tab # column 3 contains sequence identity.
awk -F'[|\t]' '{print $1","$4","$7}' temp/alignments.tab >ilocmrna_mmseqs_edgelist.csv


rm -rf temp