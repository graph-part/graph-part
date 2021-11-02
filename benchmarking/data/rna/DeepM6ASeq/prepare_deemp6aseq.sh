#!/bin/bash

# Combine individual fasta files (from deepm6aseq github) and format headers.
# https://github.com/rreybeyb/DeepM6ASeq/tree/master/data
mkdir deepm6aseq_temp
export TEMPDIR=deepm6aseq_temp
awk 'BEGIN{RS=">";OFS=";"}NR>1{print $1,$2}' dr/test_neg.fa | awk -F';' '{print ">"$2"|label=DR_NEG\n"$3}' >$TEMPDIR/dr_neg_test
awk 'BEGIN{RS=">";OFS=";"}NR>1{print $1,$2}' dr/test_pos.fa | awk -F';' '{print ">"$2"-"$3"|label=DR_POS\n"$6}'>$TEMPDIR/dr_pos_test
awk 'BEGIN{RS=">";OFS=";"}NR>1{print $1,$2}' dr/train_neg.fa | awk -F';' '{print ">"$2"|label=DR_NEG\n"$3}' >$TEMPDIR/dr_neg_train
awk 'BEGIN{RS=">";OFS=";"}NR>1{print $1,$2}' dr/train_pos.fa | awk -F';' '{print ">"$2"-"$3"|label=DR_POS\n"$6}'>$TEMPDIR/dr_pos_train

cd $TEMPDIR
cat dr_neg_test dr_pos_test dr_neg_train dr_pos_train > deepm6aseq_dr.fasta
sed 's/':'/'-'/g' deepm6aseq_dr.fasta > ../deepm6aseq_dr.fasta # fix : in headers, are only allowed as seperators.
cd -
rm -rf $TEMPDIR