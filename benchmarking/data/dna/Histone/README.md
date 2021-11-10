Original data is from https://www.cell.com/fulltext/S0092-8674(05)00645-8, but does not seem to be available directly.  
Take it from here instead:  https://github.com/Doulrs/Hilbert-CNN  


We use H3 as the task for benchmarking.
```
awk 'BEGIN{RS=">"}{if (NR!=1) {print ">"$1"|label="$3"\n"$2}}' H3.txt | sed 'y/,"/__/' >H3.fasta
```