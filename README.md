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
- Figure out how to compute sequence identity (Can get perfect match count from output, normalize ourselves? Default is full alignment length)
- Figure out needleall parameters (default BLOSUM50, penalties different from ggsearch36)
- Expose needleall parameters to CLI?

## API

Long                    | Short | Description
------------------------|-------|------------
`--fasta-file`          |`-ff`  | Path to the input fasta file, formatted according to [the input format](#Input-format).
`--out-file`            |`-of`  | Path at which to save the partition assignments as `.csv`
`--threshold`           |`-th`  | The desired partitioning threshold, should be within the bounds defined by the metric.
`--partitions`          |`-pa`  | Number of partitions to generate.
`--transformation`      |`-tf`  | Transformation to apply to the similarity/distance metric. Graph-Part operates on distances, therefore similarity metrics need to be transformed. Can be any of `one-minus`, `inverse`, `square`, `log`, `None`. See the [source](graph_part/transformations.py) for definitions. As an example, when operating with sequence identities ranging from 0 to 1, the transformation `one-minus` yields corresponding distances.
`--priority-name`       |`-pn`  | The name of the priority in the meta file. TODO what does this do
`--labels-name`         |`-ln`  | The name of the label in the meta file. Used for balancing partitions.
`--initialization-mode` |`-im`  | Use either slow or fast restricted nearest neighbor linkage or no initialization. Can be any of `slow-nn`, `fast-nn`, `simple`. Defaults to `slow-nn`.
`--threads`             |`-nt`  | The number of threads to run in parallel.
`--chunks`              |`-nc`  | The number of chunks into which to split the fasta file for multithreaded alignment.
`--load_checkpoint_path`|`-lc`  | Optional path to save the generated graph. Defaults to `None`.
`--save_checkpoint_path`|`-sc`  | Optional path to a previously generated graph for quickstart. If provided, no alignment will be performed and all arguments relating to this step are ignored.
`--edge-file`           |`-ef`  | Optional path to a comma separated file containing precomputed pairwise metrics, the first two columns should contain sequence identifiers specified in the  `--fasta-file`. This is can be used to run Graph-Part with an alignment tool different from the default `needleall`.
`--metric-column`       |`-mc`  | When using `--edge-file`, specifies in which column the metric is found. Indexing starts at 0, defaults to 2 when left unspecified.

**Flags** 
Long                    | Short | Description
------------------------|-------|------------
`--no-moving`           |`-nm`  | By default, the removing procedure tries to relocate sequences to another partition if it finds more within-threshold neighbours in any. This flag disallows moving.
`--remove-same`         |`-rs`  | This here is the inverse of removal_type (has default True), not sure what it does TODO

## FAQ
WIP
- **Can I see the process of the tool ?**  
Progress bars are tedious to implement when calling external programs, as we are doing it with `needleall`. As an alternative, when running Graph-Part in parallel mode, the progress can be inspected via `htop` and other utilies. The active `needleall` processes indicate which chunks are being processed at the moment. Chunks are processed starting from 0 up to `chunks`.

- **How should I pick `chunks` ?**  
`chunks` should be picked so that all `threads` are utilized. Each chunk is aligned to each other chunk, so `threads` <= `chunks`*`chunks` results in full utilization.

- **I want to test multiple thresholds - How can I do this efficiently ?**  
When constructing the graph, we only retain distances that are smaller than the selected `threshold`, as only those form relevant edges for partitioning the data. All other distances are discarded as they are computed. To test multiple thresholds, the most efficient way is to first try the highest threshold to be considered (when working with sequence identities, this means the lowest sequence identity) and activate checkpointing of the graph by specifiying `--save-checkpoint-path`. In the next run, use `--load-checkpoint-path` to start from your saved graph and avoid recomputing the edges. 

