'''
Adapted `hobohm2.py to have a proper script like the other baselines.

Reuse graph_part alignment code to get the edge list.
'''
from sklearn.model_selection import StratifiedKFold
import argparse
import subprocess
import os
import re
import shutil
from typing import List, Tuple
import pathlib
import json
import time
from tqdm.auto import tqdm

from graph_part.graph_part import make_graphs_from_sequences



def hobohm_homology_reduce(entity_fp: str, threshold: float = 0.3, is_nucleotide: bool = False, threads=12, alignment_mode='needle', json_dict=None)-> List[str]:


# def compute_edges(query_fp: str,
#                   library_fp: str,
#                   transformation: str,
#                   threshold: float,
#                   seq_lens: Dict[str,int],
#                   denominator = 'full',
#                   delimiter: str = '|',
#                   is_nucleotide: bool = False,
#                   gapopen: float = 10,
#                   gapextend: float = 0.5,
#                   endweight: bool = False,
#                   endopen: float = 10,
#                   endextend: float = 0.5,
#                   matrix: str = 'EBLOSUM62',
    one_minus = lambda x: 1-x
    config = {
        "alignment_mode": alignment_mode,
        "fasta_file": entity_fp,
        "threshold": one_minus(0),
        "transformation": 'one-minus',
        "priority_name": None,
        "labels_name": None,
        "denominator": 'full' if alignment_mode == 'needle' else 'shortest',
        "nucleotide": is_nucleotide,
        "prefilter": False,
        "triangular": False,
        "threads": threads,
        "chunks": 10,
        "parallel_mode": 'multiprocess',
        "gapopen": 10.0,
        "gapextend": 0.5,
        "endweight": False,
        "endopen": 10.0,
        "endextend": 0.5,
        "matrix": 'EBLOSUM62',

        # these are not used, but are required by the graph_part code.
        'partitions': 5,
    }

    # needleall -auto -stdout -aformat pair -gapopen 10.0 -gapextend 0.5 -endopen 10.0 -endextend 0.5 -datafile EBLOSUM62 -sprotein1 -sprotein2 graphpart_5.fasta.tmp graphpart_4.fasta.tmp

    if json_dict is None:
        json_dict = {}

    # get edge list
    full_graph, part_graph, labels = make_graphs_from_sequences(config, 1, json_dict, True)

    neighborlist = {}

    #debugging
    # # make an edge dataframe
    # import pandas as pd
    # edge_df = pd.DataFrame(full_graph.edges(data=True), columns=['seq1', 'seq2', 'metric'])

    # Each line in file has format: seqid_1,seqid_2,similarity
    for seq1, seq2, data in tqdm(full_graph.edges(data=True)):
        similarity = one_minus(data['metric'])

    
        # Add sequence names as we go along
        if not seq1 in neighborlist:
            neighborlist[seq1]=set()
        if not seq2 in neighborlist:
            neighborlist[seq2]=set()

        # Build lists of neighbors as we go along.
        # Note: Set.add() method automatically enforces member uniqueness - saves expensive test!
        if similarity > threshold:
            neighborlist[seq1].add(seq2)
            neighborlist[seq2].add(seq1)

    ########################################################################################

    # Build dictionary keeping track of how many neighbors each sequence has
    nr_dict = {}
    for seq in neighborlist.keys():
        nr_dict[seq]=len(neighborlist[seq])

    # Find max number of neighbors
    maxneighb = max(nr_dict.values())

    # While some sequences in list still have neighbors: remove the one with most neighbors, update counts
    # Note: could ties be dealt with intelligently?
    while maxneighb > 0:

        # Find an entry that has maxneighb neighbors, and remove it from list
        for remove_seq in nr_dict.keys():
            if nr_dict[remove_seq] == maxneighb: break
        del(nr_dict[remove_seq])

        # Update neighbor counts
        for neighbor in neighborlist[remove_seq]:
            if neighbor in nr_dict:
                nr_dict[neighbor] -= 1

        # Find new maximum number of neighbors
        maxneighb = max(nr_dict.values())

    ##############################################################################################
    # Postprocess: reinstate skipped sequences that now have no neighbors
    # Note: order may have effect. Could this be optimized?

    allseqs=set(neighborlist.keys())
    keepseqs=set(nr_dict.keys())
    skipseqs=allseqs - keepseqs

    for skipped in skipseqs:
        # if skipped sequence has no neighbors in keeplist
        if not (neighborlist[skipped] & keepseqs):
            keepseqs.add(skipped)

    return list(keepseqs)


def get_labels(identifiers: List[str], labels_name: str = 'label') -> List[str]:
    labels = []
    for ident in identifiers:

        spl = ident.strip().split('|')
        label = '0'
        for s in spl[1:]:
            if '=' in s:
                param_spl = s.split('=')
                if param_spl[0] == labels_name:
                    label = str(param_spl[1].strip())

        labels.append(label)

    return labels

def main() -> None:

    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta-file', type=str, help = 'Path to fasta file. needs Graph-Part header format.')
    parser.add_argument('--out-file', type=str, default = 'homology_reduced_split.csv')
    parser.add_argument('--labels-name', default=None)
    parser.add_argument('-th', '--threshold', type=float, default = 0.3)
    parser.add_argument('-pa', '--partitions', type=int, default=5)
    parser.add_argument('-nu', '--nucleotide', action='store_true')
    parser.add_argument('-nt', '--threads', type=int, default=12)
    parser.add_argument('-al', '--alignment-mode', type=str, default='needle')


    args = parser.parse_args()

    out_dict = {}
    out_dict['time_script_start'] = time.perf_counter()

    representatives =  hobohm_homology_reduce(args.fasta_file, args.threshold, args.nucleotide, args.threads, args.alignment_mode, out_dict)
    out_dict['time_clustering_complete'] = time.perf_counter()

    labels =  get_labels(representatives, args.labels_name)
    accs =  [x.split('|')[0] for x in representatives]

    kfold = StratifiedKFold(n_splits=args.partitions)
    folds = []

    with open(args.out_file, 'w') as f:
        f.write('AC,label-val,cluster\n')
        for i, fold in enumerate(kfold.split(representatives, labels)):
            train_idx, test_idx = fold # we only want to record the partitions, so use test_idx.
            
            for idx in test_idx:
                f.write(f'{accs[idx]},{labels[idx]},{i}\n')

    out_dict['time_script_complete'] = time.perf_counter()
    json.dump(out_dict, open(os.path.splitext(args.out_file)[0]+'_report.json','w'))



if __name__ == '__main__':
    main()