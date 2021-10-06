# Graph-Part
Protein dataset partitioning pipeline (Gíslason 2021)


# Instructions
WIP  
Parallel example command
```
python3 graphpart_ggs.py  --meta-file netgpi_dataset.fasta   --threshold 0.3 --transformation one-minus --out-file graphpart_assignments.csv --labels-name label --partitions 5 --parallel 1
```

## Input format
WIP


# TODO

- Expose multithreading n_workers to CLI - get rid of `parallel` argument, just make this n_workers>1
- Rewrite arguments processing to argparse, more concise
- Set up as python package
- Expose graph checkpointing paths/options to CLI
- Write tool to make fasta from assignments
- Implement separator choice in ggsearch_utils (hardcoded |)
- Handle > removal better, sometimes ggsearch36 parsing still includes a >. Need to remove ALL. This also implies that >at any other pos are prohibited
- Fix alignment of header in removal output:  
Min-threshold    #Entities       #Edges          Connectivity    #Problematics   #Relocated      #To-be-removed  
0.01             3539            411624                  460915                  3517            1856            1  
- Merge graphpart.py and graphpart_ggs.py . Enable previous starting from edge list as a CLI arg in the new script.

# FAQ
WIP
- **Can I see the process of the tool ?**  
Progress bars are tedious to implement when calling external programs, as we are doing it with `ggsearch36`. As an alternative, when running Graph-Part in parallel mode, the progress can be inspected via `htop` and other utilies. The active `ggsearch36` processes indicate which chunks are being processed at the moment. Chunks are processed starting from 0 up to `n_chunks`.

- **How should I pick `n_chunks` ?**  
In general, a higher number of chunks is preferable as it leads to less redundant computations (A sequence in one chunk is only aligned against each sequence in any other chunk once, but sequences within one chunk are compaired twice - align A to B and align B to A. This is a limitation as `ggsearch36` was not designed to compute distance matrices).

- **I want to test multiple thresholds - How can I do this efficiently ?**  
When constructing the graph, we only retain distances that are smaller than the selected `threshold`, as only those form relevant edges for partitioning the data. All other distances are discarded as they are computed. To test multiple thresholds, the most efficient way is to first try the highest threshold to be considered (when working with sequence identities, this means the lowest sequence identity) and activate checkpointing of the graph.
