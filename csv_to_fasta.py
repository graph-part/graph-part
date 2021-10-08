'''
Utility to convert .csv datasets into
Graph-Part compatible fasta files.

TODO
allow more than 1 column for id, label, prio? 
figure out how to best include in package/cli
'''

import pandas as pd
import argparse



def convert_to_fasta(csv_file: str, 
                    sequence_col:int, 
                    identifier_col: int, 
                    label_col: int = None, 
                    priority_col: int = None, 
                    label_name: str = 'label',
                    priority_name: str='priority',
                    delimiter='|'):

    df = pd.read_csv(csv_file)

    cols = df.columns

    seqs = df[cols[sequence_col]]

    # make fasta headers
    headers = '>' + df[cols[identifier_col]]
    if label_col is not None:
        headers = headers + delimiter + label_name + '=' + df[cols[label_col]]
    if priority_col is not None:
        headers = headers + delimiter + priority_name + '=' + df[cols[priority_col]]

    target_path = csv_file[:-3] + 'fasta'
    with open(target_path, 'w') as f:
            for head, seq in zip(headers, seqs):
                f.write(head+'\n')
                f.write(seq+'\n')

    print(f'Created {target_path}')



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', '-f', type=str, help='Path the the .csv file to be converted.')
    parser.add_argument('--sequence_col', '-s', type=int, help='Index of the column that contains the amino acid sequence.')
    parser.add_argument('--identifier_col', '-id', type=int, help='Index of the column that contains the sequence identifier')
    parser.add_argument('--label_col', '-la', type=int, default=None, help='Index of the column that contains the label of the sequence.')
    parser.add_argument('--priority_col', '-pr', type=int, default=None, help='Index of the column that contains the priority group of the sequence.')
    parser.add_argument('--delimiter', '-d', type=str, default='|', help='Delimiter to use in the generated fasta headers.')
    parser.add_argument('--label_name', '-ln', type=str, default='label', help='Name of the label in the generated headers.')
    parser.add_argument('--priority_name', '-pn', type=str, default='priority', help='Name of the priority in the generated headers.')
    args = parser.parse_args()

    convert_to_fasta(args.file, 
                     args.sequence_col, 
                     args.identifier_col, 
                     args.label_col, 
                     args.priority_col, 
                     args.label_name, 
                     args.priority_name, 
                     args.delimiter
                     )

if __name__ == '__main__':
    main()