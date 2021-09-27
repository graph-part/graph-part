'''
Helper script to generate a fasta file from .csv so 
that it can be used by the partitioning pipeline.
'''

import pandas as pd
import argparse
import os


def csv_to_fasta(in_path, out_path, identifier_col, label_col, sequence_col):
    
    df = pd.read_csv(in_path)

    with open(out_path, 'w') as f:
        for idx, row in df.iterrows():
            f.write('>'+row[identifier_col].replace(' ', '')+'|label='+row[label_col]+'\n')
            f.write(row[sequence_col][:70]+'\n')



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out_path', type=str, help='path to output file')
    parser.add_argument('--in_path', type=str, help='path to input file')
    parser.add_argument('--identifier_col', type=str, help='Sequence ID col in csv. Whitespace will be removed.', default='Entry')
    parser.add_argument('--sequence_col', type=str, help='AA sequence col in csv', default='Sequence')
    parser.add_argument('--label_col', type=str, help='Class label col in csv')

    args = parser.parse_args()

    csv_to_fasta(args.in_path, args.out_path, args.identifier_col, args.label_col, args.sequence_col)