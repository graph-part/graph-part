# Graph-Part
Protein dataset partitioning pipeline (Gíslason 2021)

## Installation

Graph-Part relies on [needleall](https://www.bioinformatics.nl/cgi-bin/emboss/help/needleall) from the [EMBOSS](http://emboss.sourceforge.net/) package for fast Needleman-Wunsch alignments of sequences. Please refer to the official EMBOSS documentation for installation methods.

We recommend to install Graph-Part in a conda environment, and install EMBOSS from [bioconda](https://anaconda.org/bioconda/emboss) via 
```
conda install -c bioconda emboss
```

Alternatively, on Ubuntu, EMBOSS is available directly via `sudo apt-get install emboss` .

###TODO pip graphpart

## Instructions
WIP  
Parallel example command
```
graphpart  --fasta-file netgpi_dataset.fasta --threshold 0.3 --transformation one-minus --out-file graphpart_assignments.csv --labels-name label --partitions 5 --threads 12
```

## Input format
WIP


# TODO

- Write tool to make fasta from assignments (concatenate assignment to header, seperator of choice)
- Write tool to make fasta from .csv (specifify spearator, label_col (multiple?) and priority_col)
- Fix alignment of header in removal output - this seems to happen with large Connectivity values:  
```
Min-threshold    #Entities       #Edges          Connectivity    #Problematics   #Relocated      #To-be-removed  
0.01             3539            411624                  460915                  3517            1856            1  
```


## FAQ
WIP
- **Can I see the process of the tool ?**  
Progress bars are tedious to implement when calling external programs, as we are doing it with `needleall`. As an alternative, when running Graph-Part in parallel mode, the progress can be inspected via `htop` and other utilies. The active `needleall` processes indicate which chunks are being processed at the moment. Chunks are processed starting from 0 up to `chunks`.

- **How should I pick `chunks` ?**  
`chunks` should be picked so that all `threads` are utilized. Each chunk is aligned to each other chunk, so `threads` <= `chunks`*`chunks` results in full utilization.

- **I want to test multiple thresholds - How can I do this efficiently ?**  
When constructing the graph, we only retain distances that are smaller than the selected `threshold`, as only those form relevant edges for partitioning the data. All other distances are discarded as they are computed. To test multiple thresholds, the most efficient way is to first try the highest threshold to be considered (when working with sequence identities, this means the lowest sequence identity) and activate checkpointing of the graph by specifiying `--save-checkpoint-path`. In the next run, use `--load-checkpoint-path` to start from your saved graph and avoid recomputing the edges. 

