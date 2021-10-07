'''
This script calls ggsearch36 to compute pairwise identities between all sequences and parses the outputs. 
Generates and edge list with pairwise % identity.
'''

import subprocess
from os import path, makedirs
import argparse

# edgelist_out = '/work3/felteu/edgelist_rerun_14122020.csv'
# #when self-aligning huge datasets, use multiple batches. query is subset of data, lib is full data.
# #cat the files after processing, graph_part does not care about order when parsing.
# fasta_file_query = 'signalp_6_seqs_only_graphpartheaders.fasta'#'/work3/felteu/edgelist_splitted/signalp6_seqonly_for_graphpart.fasta'#full_updated_data_seqs_only.fasta'
# fasta_file_lib = 'signalp_6_seqs_only_graphpartheaders.fasta'#'/work3/felteu/edgelist_splitted/signalp6_seqonly_for_graphpart.fasta'
# ggs = path.expanduser('/work3/felteu/fasta36/bin/ggsearch36')

def run_ggsearch(fasta_file, ggsearch_path, output_file, ref_file):
    ggs = path.expanduser(ggsearch_path)
    with subprocess.Popen(
            [ggs,"-E","41762",fasta_file,ref_file],
            stdout=subprocess.PIPE,
            bufsize=1,
            universal_newlines=True) as proc:
        with open(output_file,'w+') as outf:
            for line_nr, line in enumerate(proc.stdout):
                if '>>>' in line:
                    qry_nr = int(line[2])
                    this_qry = line[6:70].split()[0].split('|')[0]

                elif line[0:2] == '>>':
                    this_lib = line[2:66].split()[0].split('|')[0]

                elif line[:13] == 'global/global':
                    identity = float(line.split()[4][:-1])/100
                    #print(qry_nr, this_qry, this_lib, identity)

                    if this_qry == this_lib:
                        continue
                    outf.write("%s,%s,%.3f\n" % (this_qry, this_lib,identity))



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--sequences', type=str, help='.fasta file of sequences to be compared.')
    parser.add_argument('--sequences_2', type=str, help='.fasta file of sequences to be compared.')
    parser.add_argument('--outfile', type=str, help='path of output file to be created.')
    parser.add_argument('--ggs', type=str, help='Path to ggsearch36 executable.')

    args = parser.parse_args()

    run_ggsearch(args.sequences, args.ggs, args.outfile, args.sequences_2)