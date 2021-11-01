Download here https://services.healthtech.dtu.dk/service.php?SignalP-5.0  


```
# remove 3rd lines and prefix label with "label="
awk '{if (NR%3!=0){print $0}}' train_set.fasta | sed 's/|/|label=/2' >train_set_reformatted.fasta
```