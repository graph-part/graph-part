from sklearn.model_selection import StratifiedKFold
import argparse
import subprocess
import os
import re
import shutil
from typing import List, Tuple
import time
import json


def mmseqs2_homology_reduce(entity_fp: str, threshold: float = 0.3, cluster_mode: int = 0, seq_id_mode: int = 0)-> List[str]:

    out_prefix = 'reduction_result'
    temp_path = 'mmseqs_temp'
    os.makedirs(temp_path, exist_ok=True)

    subprocess.run(['mmseqs', 'easy-cluster', '--min-seq-id', str(threshold), '--seq-id-mode', str(seq_id_mode), '--cluster-mode', str(cluster_mode), entity_fp, out_prefix, temp_path])
    shutil.rmtree(temp_path)

    representatives = []
    with open(f'{out_prefix}_cluster.tsv') as f:
        for line_nr, line in enumerate(f):
            spl = line.strip().split('\t')

            rep = spl[0]
            representatives.append(rep)

    for x in os.listdir():
        if 'reduction_result' in x:
            os.remove(x)

    return list(set(representatives))


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

    parser.add_argument('--cluster-mode', type=int, default=0)
    parser.add_argument('--seq-id-mode', type=int, default=0)

    args = parser.parse_args()

    out_dict = {}
    out_dict['time_script_start'] = time.perf_counter()

    representatives =  mmseqs2_homology_reduce(args.fasta_file, args.threshold, args.cluster_mode, args.seq_id_mode)
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