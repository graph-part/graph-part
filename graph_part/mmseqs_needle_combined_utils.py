'''
This file contains functions that first call mmseqs2 to compute the 
pairwise sequence identities and insert them into the graph.
Identities that are close to the trehshold are then recomputed with
needleall to get the exact metric value.

This is faster than running needleall on the full dataset.
'''
import networkx as nx
import os
import tempfile
import numpy as np
from typing import Dict, List, Tuple
import concurrent.futures
import subprocess
from tqdm.auto import tqdm
from .transformations import TRANSFORMATIONS
from .mmseqs_utils import generate_edges_mmseqs
from .needle_utils import parse_fasta, NORMALIZATIONS


# the elementary alignment operation for a pair. 
def _align_two_single_sequence_files(
    file_1: str,
    file_2: str,
    transformation: str,
    seq_lens: Dict[str, int],
    denominator: str = 'full',
    is_nucleotide: bool = False,
    gapopen: float = 10,
    gapextend: float = 0.5,
    endweight: bool = True,
    endopen: float = 10,
    endextend: float = 0.5,
    matrix: str = 'EBLOSUM62',
    ) -> float:
    '''Runs needle alignment on two single sequence files and returns the transformed identity.'''

    if is_nucleotide:
        type_1, type_2, = '-snucleotide1', '-snucleotide2'
    else:
        type_1, type_2 = '-sprotein1', '-sprotein2'

    command = ["needleall","-auto","-stdout", 
            "-aformat", "pair", 
            "-gapopen", str(gapopen),
            "-gapextend", str(gapextend),
            "-endopen", str(endopen),
            "-endextend", str(endextend),
            "-datafile", matrix,
            type_1, type_2, file_1, file_2]
    if endweight:
        command = command + ["-endweight"]

    count = 0
    import subprocess
    with subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            bufsize=1,
            universal_newlines=True) as proc:
        for line_nr, line in enumerate(proc.stdout):  
            if  line.startswith('# 1:'):

                # # 1: P0CV73
                this_qry = line[5:].split()[0].split('|')[0]

            elif line.startswith('# 2:'):
                this_lib = line[5:].split()[0].split('|')[0]

            elif line.startswith('# Identity:'):
                identity_line = line

            elif line.startswith('# Gaps:'):
                count +=1

                gaps, rest = line[7:].split('/')
                gaps = int(gaps)
                length = int(rest.split('(')[0])

                
                # Compute different sequence identities as needed.
                if denominator == 'full': # full is returned by default. just need to parse
                    identity = float(identity_line.split('(')[1][:-3])/100
                elif denominator == 'no_gaps':
                    n_matches =  int(identity_line[11:].split('/')[0]) #int() does not mind leading spaces
                    identity = float(n_matches/(length-gaps))
                else:
                    n_matches =  int(identity_line[11:].split('/')[0]) #int() does not mind leading spaces
                    identity = NORMALIZATIONS[denominator](n_matches, seq_lens[this_qry], seq_lens[this_lib])
                
                try:
                    metric = TRANSFORMATIONS[transformation](identity)
                except ValueError or TypeError:
                    raise TypeError("Failed to interpret the identity value %r. Please ensure that the ggsearch36 output is correctly formatted." % (identity))
                
                return metric
    
    # if we reach this point, needleall did not return an alignment
    # This happens if the minimum score is not reached, this
    # will be reported in needleall.error (even if it isn't really an error)
    return TRANSFORMATIONS[transformation](0.0)


# this handles lists of pairs
# def _generate_needle_edges_pairwise(
#         sequences: Dict[str,str],
#         pairs: List[Tuple[str, str]],
#         transformation: str,
#         denominator: str = 'full',
#         triangular: bool = False,
#         is_nucleotide: bool = False,
#         gapopen: float = 10,
#         gapextend: float = 0.5,
#         endweight: bool = True,
#         endopen: float = 10,
#         endextend: float = 0.5,
#         matrix: str = 'EBLOSUM62',
#         ) -> Dict[Tuple[str,str], float]:
#     '''
#     We normally call needleall to compute all-vs-all sequence identities on whole
#     files. Now, we only need to calculate pre-selected pairs. This function
#     takes a list of pairs and computes the sequence identities for each pair
#     using needleall. The results are then inserted into the graph.
#     '''
#     out_dict = {}
#     # create a temporary directory
#     with tempfile.TemporaryDirectory() as temp_dir:
#         for n1, n2 in pairs:
#             seq1 = sequences[n1]
#             seq2 = sequences[n2]
#             # write the sequences to a temporary file
#             with open(os.path.join(temp_dir, f'{n1}.fasta'), 'w') as f:
#                 f.write(f'>{n1}\n{seq1}\n')
#             with open(os.path.join(temp_dir, f'{n2}.fasta'), 'w') as f:
#                 f.write(f'>{n2}\n{seq2}\n')
#             # run needleall
#             metric = _align_two_single_sequence_files(
#                 file_1 = os.path.join(temp_dir, f'{n1}.fasta'),
#                 file_2 = os.path.join(temp_dir, f'{n2}.fasta'),
#                 transformation = transformation,
#                 seq_lens = {n1: len(seq1), n2: len(seq2)},
#                 denominator = denominator,
#                 is_nucleotide = is_nucleotide,
#                 gapopen = gapopen,
#                 gapextend = gapextend,
#                 endweight = endweight,
#                 endopen = endopen,
#                 endextend = endextend,
#                 matrix = matrix,
#             )
#             if not triangular:
#                 # also run the other direction
#                 metric_2 = _align_two_single_sequence_files(
#                     file_1 = os.path.join(temp_dir, f'{n2}.fasta'),
#                     file_2 = os.path.join(temp_dir, f'{n1}.fasta'),
#                     transformation = transformation,
#                     seq_lens = {n1: len(seq1), n2: len(seq2)},
#                     denominator = denominator,
#                     is_nucleotide = is_nucleotide,
#                     gapopen = gapopen,
#                     gapextend = gapextend,
#                     endweight = endweight,
#                     endopen = endopen,
#                     endextend = endextend,
#                     matrix = matrix,
#                 )
#                 metric = min(metric, metric_2)

#             out_dict[(n1, n2)] = metric

#     return out_dict


# an alternative implementation of the above function that does not make single sequence files, but just
# runs a batch of them directly through needleall.
# this computes alignments we don't need, but can be faster if the overhead of creating the files is large.
def _generate_needle_edges_pairwise_batch(
        qry_file: str,
        lib_file: str,
        # sequences: Dict[str,str],
        # pairs: List[Tuple[str, str]],
        seq_lens: Dict[str, int],
        transformation: str,
        denominator: str = 'full',
        triangular: bool = False,
        is_nucleotide: bool = False,
        gapopen: float = 10,
        gapextend: float = 0.5,
        endweight: bool = True,
        endopen: float = 10,
        endextend: float = 0.5,
        matrix: str = 'EBLOSUM62',
        ) -> Dict[Tuple[str,str], float]:
    '''
    We normally call needleall to compute all-vs-all sequence identities on whole
    files. Now, we only need to calculate pre-selected pairs. This function
    takes a list of pairs and computes the sequence identities for each pair
    using needleall. The results are then inserted into the graph.
    '''
    out_dict = {}


    # run needleall
    if is_nucleotide:
        type_1, type_2, = '-snucleotide1', '-snucleotide2'
    else:
        type_1, type_2 = '-sprotein1', '-sprotein2'
    command = ["needleall","-auto","-stdout",
            "-aformat", "pair",
            "-gapopen", str(gapopen),
            "-gapextend", str(gapextend),
            "-endopen", str(endopen),
            "-endextend", str(endextend),
            "-datafile", matrix,
            type_1, type_2,
            qry_file,
            lib_file]
    if endweight:
        command = command + ["-endweight"]

    count = 0
    with subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        bufsize=1,
        universal_newlines=True) as proc:

        for line_nr, line in enumerate(proc.stdout):  
            if  line.startswith('# 1:'):

                # # 1: P0CV73
                this_qry = line[5:].split()[0].split('|')[0]

            elif line.startswith('# 2:'):
                this_lib = line[5:].split()[0].split('|')[0]

            elif line.startswith('# Identity:'):
                identity_line = line

            elif line.startswith('# Gaps:'):
                count +=1

                gaps, rest = line[7:].split('/')
                gaps = int(gaps)
                length = int(rest.split('(')[0])

                
                # Compute different sequence identities as needed.
                if denominator == 'full': # full is returned by default. just need to parse
                    identity = float(identity_line.split('(')[1][:-3])/100
                elif denominator == 'no_gaps':
                    n_matches =  int(identity_line[11:].split('/')[0]) #int() does not mind leading spaces
                    identity = float(n_matches/(length-gaps))
                else:
                    n_matches =  int(identity_line[11:].split('/')[0]) #int() does not mind leading spaces
                    identity = NORMALIZATIONS[denominator](n_matches, seq_lens[this_qry], seq_lens[this_lib])
                
                try:
                    metric = TRANSFORMATIONS[transformation](identity)
                except ValueError or TypeError:
                    raise TypeError("Failed to interpret the identity value %r. Please ensure that the ggsearch36 output is correctly formatted." % (identity))
                
                out_dict[(this_qry, this_lib)] = metric

    return out_dict
        
                                

def generate_needle_edges_pairwise_mp(
        sequences: Dict[str,str],
        pairs: List[Tuple[str, str]],
        full_graph: nx.classes.graph.Graph, 
        transformation: str,
        threshold: float,
        denominator: str = 'full',
        n_procs: int = 4,
        parallel_mode: str = 'multithread',
        triangular: bool = False,
        delimiter: str = '|',
        is_nucleotide: bool = False,
        gapopen: float = 10,
        gapextend: float = 0.5,
        endweight: bool = True,
        endopen: float = 10,
        endextend: float = 0.5,
        matrix: str = 'EBLOSUM62', 
        ) -> None:
    '''
    Splits up the list of pairs so that each parallel worker handles one.
    Collects the results and inserts them into the graph.
    '''

    # split the list of pairs
    chunks = np.array_split(pairs, len(pairs)/1)

    # create a temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        for idx, pairs in tqdm(enumerate(chunks), desc='Writing files', leave=False):

            # write the sequences to files for alignment
            query_file = open(os.path.join(temp_dir, f'{idx}_query.fasta'), 'w')
            lib_file = open(os.path.join(temp_dir, f'{idx}_lib.fasta'), 'w')
            for n1, n2 in pairs:
                query_file.write(f'>{n1}\n{sequences[n1]}\n')
                lib_file.write(f'>{n2}\n{sequences[n2]}\n')

            query_file.close()
            lib_file.close()

        # prepare seq lens once for all workers
        seq_lens = {x: len(y) for x, y in sequences.items()}

        # create a pool of workers
        if parallel_mode == 'multithread':
            executor = concurrent.futures.ThreadPoolExecutor(n_procs)
        elif parallel_mode == 'multiprocess':
            executor = concurrent.futures.ProcessPoolExecutor(n_procs)
        else:
            raise ValueError(f'Unknown parallel mode {parallel_mode}')

        # print(f'Each chunk is {len(chunks[0])} pairs.')
        
        # submit the jobs
        futures = []
        for idx, chunk in enumerate(chunks):
            futures.append(executor.submit(
                _generate_needle_edges_pairwise_batch,#_generate_needle_edges_pairwise,
                qry_file = os.path.join(temp_dir, f'{idx}_query.fasta'),
                lib_file = os.path.join(temp_dir, f'{idx}_lib.fasta'),
                seq_lens = seq_lens,
                transformation = transformation,
                denominator = denominator,
                triangular = triangular,
                is_nucleotide = is_nucleotide,
                gapopen = gapopen,
                gapextend = gapextend,
                endweight = endweight,
                endopen = endopen,
                endextend = endextend,
                matrix = matrix,
            ))

        # collect the results
        pbar = tqdm(total=len(pairs))
        for future in concurrent.futures.as_completed(futures):
            if future.exception() is not None:
                print(future.exception())
                raise RuntimeError('One of the alignment processes did not complete sucessfully.')
            out_dict = future.result()
            for pair, metric in out_dict.items():
                
                # if the exact metric does not violate the threshold,
                # we can remove it from the graph
                # we need to check whether it's actually in the graph though,
                # as in chunk mode we also compute some alignments that are not
                # required. needeall limitations ...
                if (metric > threshold) and full_graph.has_edge(pair[0], pair[1]):
                    full_graph.remove_edge(pair[0], pair[1])

                # otherwise, we insert the metric into the graph
                else:
                    full_graph.add_edge(pair[0], pair[1], metric=metric, delimiter=delimiter)
                    

            pbar.update(len(out_dict))

        # shutdown the pool of workers
        executor.shutdown(wait=True)


def generate_edges_mmseqs_needle_combined(
        entity_fp: str, 
        full_graph: nx.classes.graph.Graph, 
        transformation: str,
        threshold: float,
        recompute_threshold: float,
        denominator_needle: str = 'full',
        denominator_mmseqs: str = 'longest',
        n_procs: int = 4,
        parallel_mode: str = 'multithread',
        triangular: bool = False,
        delimiter: str = '|',
        is_nucleotide: bool = False,
        gapopen: float = 10,
        gapextend: float = 0.5,
        endweight: bool = True,
        endopen: float = 10,
        endextend: float = 0.5,
        matrix: str = 'EBLOSUM62',
        use_prefilter: bool = False,
        ) -> None:
    '''
    First we run mmseqs2 on all sequences.
    If their mmseqs2 identity is below the recompute_threshold,
    we recompute the exact value with needleall.

    The recompute_threshold therefore needs to be higher than the
    final threshold.
    '''
    
    # first run mmseqs2 to get the approximate values
    # in this case, all values need to be inserted into the graph
    # normally, we only keep the ones below the threshold
    generate_edges_mmseqs(
        entity_fp = entity_fp,
        full_graph = full_graph,
        tranformation = transformation,
        threshold_transformed = recompute_threshold, # NOTE: this is a dummy value, since we want to insert all values
        threshold_original = None,
        denominator = denominator_mmseqs,
        delimiter = delimiter,
        is_nucleotide = is_nucleotide,
        use_prefilter = use_prefilter,
    )

    # above command inserted all pairwise distances into the graph.
    # MMseqs distances are overestimated. We now recompute the exact values
    # for all pairs with a distance above the recompute_threshold.
    pairs = []
    for u, v, data in full_graph.edges(data=True):
        if data['metric'] > threshold:
            pairs.append((u,v))

    # get the sequences, make dict
    ids, seqs = parse_fasta(entity_fp)
    ids = [i[1:] for i in ids] # get clean ids ('>')
    seqs = dict(zip(ids, seqs))

    print('Computing NW values for %i pairs.' % len(pairs))

    # compute the exact values
    generate_needle_edges_pairwise_mp(
        sequences = seqs,
        pairs = pairs,
        full_graph = full_graph,
        transformation = transformation,
        threshold = threshold,
        denominator = denominator_needle,
        n_procs = n_procs,
        parallel_mode = parallel_mode,
        triangular = triangular,
        delimiter = delimiter,
        is_nucleotide = is_nucleotide,
        gapopen = gapopen,
        gapextend = gapextend,
        endweight = endweight,
        endopen = endopen,
        endextend = endextend,
        matrix = matrix,
    )







