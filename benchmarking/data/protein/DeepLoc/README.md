Download here https://services.healthtech.dtu.dk/service.php?DeepLoc-1.0






```
awk 'BEGIN{RS=">"}{print $1"\t"$2"\t"$3"\t"$4;}' deeploc_data.fasta >deeploc_data.tsv

grep -v "test" deeploc_data.tsv | awk '{print ">"$1"|label="$2"\n"$3}' >deeploc_train.fasta
grep "test" deeploc_data.tsv | awk '{print ">"$1"|label="$2"\n"$4}' >deeploc_test.fasta

sed -i 's/-[SMU]//' deeploc_train.fasta
sed -i 's/-[SMU]//' deeploc_test.fasta

cat deeploc_train.fasta deeploc_test.fasta >deeploc_all.fasta
```