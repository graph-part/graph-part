'''
A naive baseline that ignores homology altogether and randomly splits the data into partitions.
'''
from sklearn.model_selection import StratifiedKFold
import argparse
import os
from typing import List, Tuple
import pathlib

def parse_fasta(fasta_fp: str) -> Tuple[List[str], List[str]]:
    '''
    Parses a two-line fasta file and returns a tuple of identifiers and sequences.
    '''
    identifiers = []
    sequences = []
    with open(fasta_fp) as f:
        for line in f:
            if line.startswith('>'):
                identifiers.append(line.strip().lstrip('>'))
            else:
                sequences.append(line.strip())

    return identifiers, sequences

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

    args = parser.parse_args()

    identifiers, sequences = parse_fasta(args.fasta_file)

    labels =  get_labels(identifiers, args.labels_name)
    accs =  [x.split('|')[0] for x in identifiers]

    kfold = StratifiedKFold(n_splits=args.partitions)
    folds = []

    with open(args.out_file, 'w') as f:
        f.write('AC,label-val,cluster\n')
        for i, fold in enumerate(kfold.split(identifiers, labels)):
            train_idx, test_idx = fold # we only want to record the partitions, so use test_idx.
            
            for idx in test_idx:
                f.write(f'{accs[idx]},{labels[idx]},{i}\n')



if __name__ == '__main__':
    main()