from sklearn.model_selection import StratifiedKFold
import argparse
import subprocess
import os
import re
import shutil
from typing import List, Tuple
import pathlib




def psicdhit_homology_reduce(entity_fp: str, threshold: float = 0.3)-> List[str]:

    out_prefix = 'reduction_result'
   
    subprocess.run(['cd-hit', '-i', entity_fp, '-o', 'pre_reduced_1' , '-c', '0.9', '-n', '5', '-T', '0' ])
    subprocess.run(['cd-hit', '-i', 'pre_reduced_1', '-o', 'pre_reduced_2' , '-c', '0.6', '-n', '4', '-T', '0' ]) #1792

    psi_path = pathlib.Path(__file__).parent.resolve() / 'psi_cd_hit.pl'
        #subprocess.run(['perl','-I', psi_path.parent, psi_path,
    subprocess.run(['perl', psi_path, '-i', 'pre_reduced_2', '-o', out_prefix, '-c', str(threshold)])


    representatives = []
    with open('reduction_result') as f:
        for line in f:
            if line.startswith('>'):
                rep = line.strip().lstrip('>')
                representatives.append(rep)

    for x in os.listdir():
        if 'reduction_result' in x:
            if os.path.isdir(x):
                shutil.rmtree(x)
            else:
                os.remove(x)
        elif 'pre_reduced' in x:
            if os.path.isdir(x):
                shutil.rmtree(x)
            else:
                os.remove(x)

    return representatives



def cdhit_homology_reduce(entity_fp: str, threshold: float = 0.3, is_nucleotide: bool = False)-> List[str]:

    if not is_nucleotide and threshold<0.4:
        return psicdhit_homology_reduce(entity_fp, threshold)
    if is_nucleotide and threshold<0.4:
        raise NotImplementedError('Threshold too low for this data type.')

    out_prefix = 'reduction_result'

    word_size = '2'
    if threshold >= 0.5:
        word_size = '3'
    if threshold >= 0.6:
        word_size = '4'
    if threshold >= 0.7:
        word_size = '5'

    if is_nucleotide:
        subprocess.run(['cd-hit-est', '-i', entity_fp, '-o', out_prefix , '-c', str(threshold), '-n', word_size, '-T', '0' ])
    else:
        subprocess.run(['cd-hit', '-i', entity_fp, '-o', out_prefix , '-c', str(threshold), '-n', word_size, '-T', '0' ])

    representatives = []
    with open('reduction_result') as f:
        for line in f:
            if line.startswith('>'):
                rep = line.strip().lstrip('>')
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

    args = parser.parse_args()

    representatives =  cdhit_homology_reduce(args.fasta_file, args.threshold)

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



if __name__ == '__main__':
    main()