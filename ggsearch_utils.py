'''
This file contains functions that call ggsearch36 on a provided
.fasta file and insert the computed pairwise sequence identities
as edges into a provided networkx graph.
'''
import multiprocessing
import networkx as nx
from os import path, remove
import numpy as np
import math
from itertools import groupby
from typing import Dict, List, Tuple, Iterator
import concurrent.futures

## Metric transformations should be defined here
TRANSFORMATIONS = {
    'one-minus': lambda x: 1-x, 
    'inverse': lambda x: 1/x if x > 0 else float('Inf'), 
    'square': lambda x: x**2,
    'log': lambda x: np.log(x),
    'none': lambda x: x,
    'None': lambda x: x,
    None: lambda x: x
}

def parse_fasta(fastafile: str) -> Tuple[List[str],List[str]]:
	'''
    Parses fasta file into lists of identifiers and sequences.
	Can handle multi-line sequences and empty lines.
    '''
	ids = []
	seqs = []
	with open(fastafile, 'r') as f:

		id_seq_groups = (group for group in groupby(f, lambda line: line.startswith(">")))

		for is_id, id_iter in id_seq_groups:
			if is_id: # Only needed to find first id line, always True thereafter
			    ids.append(next(id_iter).strip())
			    seqs.append("".join(seq.strip() for seq in next(id_seq_groups)[1]))
        
	return ids, seqs



def generate_edges(entity_fp: str, 
                  full_graph: nx.classes.graph.Graph, 
                  tranformation: str,
                  threshold: float,
                  ggsearch_path: str,
                  separator: str = '|',
                  ) -> None:
    '''
    Call ggsearch36 and insert found edges into the graph as they are computed.
    This is the default implementation that runs one single process for the full
    dataset without multithreading.
    '''

    # get length of dataset as e-value for ggsearch36
    e_value = str(len(parse_fasta(entity_fp)[0]) +1 ) #string conversion required by subprocess API
    import subprocess
    ggs = path.expanduser(ggsearch_path)
    with subprocess.Popen(
            [ggs,"-E", e_value ,entity_fp, entity_fp],
            stdout=subprocess.PIPE,
            bufsize=1,
            universal_newlines=True) as proc:
        #with open(output_file,'w+') as outf:
        for line_nr, line in enumerate(proc.stdout):
            if '>>>' in line:
                qry_nr = int(line[2])
                this_qry = line[6:70].split()[0].split(separator)[0]
                this_qry = this_qry.lstrip('>')

            elif line[0:2] == '>>':
                this_lib = line[2:66].split()[0].split(separator)[0]
                this_lib = this_lib.lstrip('>')

            elif line[:13] == 'global/global':
                identity = float(line.split()[4][:-1])/100

                try:
                    metric = TRANSFORMATIONS[tranformation](identity)
                except ValueError or TypeError:
                    raise TypeError("Failed to interpret the metric column value %r. Please ensure that the edge list file is correctly formatted and that the correct column is specified." % (identity))
                
                if this_qry == this_lib:
                    continue
                if metric > threshold:
                    continue
                # NOTE this case should raise an error - graph was constructed from same file before, and so all the nodes should be there.
                if not full_graph.has_node(this_qry) or not full_graph.has_node(this_lib):
                    raise RuntimeError(f'Tried to insert edge {this_qry}-{this_lib} into the graph, but did not find nodes. This should not happen, please report a bug.')
                if full_graph.has_edge(this_qry, this_lib):
                    if full_graph[this_qry][this_lib]['metric'] > metric:
                        nx.set_edge_attributes(full_graph,{(this_qry,this_lib):metric}, 'metric')
                else:
                    full_graph.add_edge(this_qry, this_lib, metric=metric)  




def chunk_fasta_file(ids: List[str], seqs: List[str], n_chunks: int) -> int:
    '''
    Break up fasta file into multiple smaller files that can be
    used for multiprocessing.
    Returns the number of generated chunks.
    '''

    chunk_size = math.ceil(len(ids)/n_chunks)

    empty_chunks = 0
    for i in range(n_chunks):
        # because of ceil() we sometimes make less partitions than specified.
        if i*chunk_size>len(ids):
            empty_chunks +=1
            continue

        chunk_ids = ids[i*chunk_size:(i+1)*chunk_size]
        chunk_seqs = seqs[i*chunk_size:(i+1)*chunk_size]

        with open(f'graphpart_{i}.fasta.tmp', 'w') as f:
            for id, seq in zip(chunk_ids, chunk_seqs):
                f.write(id+'\n')
                f.write(seq+'\n')

    return n_chunks - empty_chunks


def compute_edges(query_fp: str,
                  library_fp: str,
                  full_graph: nx.classes.graph.Graph, 
                  transformation: str,
                  threshold: float,
                  e_value: float,
                  ggsearch_path: str,
                  separator: str = '|',
                  ) -> None:
    '''
    Run ggsearch36 on query_fp and library_fp,
    Retrieve pairwise similiarities, transform and
    insert into edge_dict.
    '''
    import subprocess
    ggs = path.expanduser(ggsearch_path)
    with subprocess.Popen(
            [ggs,"-E", e_value, query_fp, library_fp],
            stdout=subprocess.PIPE,
            bufsize=1,
            universal_newlines=True) as proc:
        for line_nr, line in enumerate(proc.stdout):
            if '>>>' in line:
                qry_nr = int(line[2])
                this_qry = line[6:70].split()[0].split(separator)[0]

            elif line[0:2] == '>>':
                this_lib = line[2:66].split()[0].split(separator)[0]

            elif line[:13] == 'global/global':
                identity = float(line.split()[4][:-1])/100
                #print(qry_nr, this_qry, this_lib, identity)

                try:
                    metric = TRANSFORMATIONS[transformation](identity)
                except ValueError or TypeError:
                    raise TypeError("Failed to interpret the identity value %r. Please ensure that the ggsearch36 output is correctly formatted." % (identity))
                
                if this_qry == this_lib:
                    continue
                if metric > threshold:
                    continue
                # NOTE this case should raise an error - graph was constructed from same file before, and so all the nodes should be there.
                if not full_graph.has_node(this_qry) or not full_graph.has_node(this_lib):
                    raise RuntimeError(f'Tried to insert edge {this_qry}-{this_lib} into the graph, but did not find nodes. This should not happen, please report a bug.')
                if full_graph.has_edge(this_qry, this_lib):
                    if full_graph[this_qry][this_lib]['metric'] > metric:
                        nx.set_edge_attributes(full_graph,{(this_qry,this_lib):metric}, 'metric')
                else:
                    full_graph.add_edge(this_qry, this_lib, metric=metric)  



def generate_edges_mp(entity_fp: str, 
                  full_graph: nx.classes.graph.Graph, 
                  transformation: str,
                  threshold: float,
                  ggsearch_path: str,
                  n_chunks: int = 10,
                  n_procs: int = 4) -> None:
    '''
    Call ggsearch36 to compute all pairwise sequence identities in the dataset.
    Uses chunked fasta files and multiple threads with ggsearch36 subprocesses 
    to speed up computation.
    '''

    # chunk the input
    ids, seqs = parse_fasta(entity_fp)
    e_value = str(len(ids)+1)

    n_chunks = chunk_fasta_file(ids, seqs, n_chunks) #get the actual number of generated chunks.

    # start n_procs threads, each thread starts a subprocess
    # Because of threading's GIL we can write edges directly to the full_graph object.
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_procs) as executor:
        for i in range(n_chunks):
            for j in range(i,n_chunks):
                # need to call each chunk with a) itself and b) all others.
                # distance is symmetric, so after having 0:1 don't need 1:0
                q = f'graphpart_{i}.fasta.tmp'
                l = f'graphpart_{j}.fasta.tmp'
                future = executor.submit(compute_edges, q, l, full_graph, transformation, threshold, e_value, ggsearch_path)

    #delete the chunks
    for i in range(n_chunks):
        remove(f'graphpart_{i}.fasta.tmp')