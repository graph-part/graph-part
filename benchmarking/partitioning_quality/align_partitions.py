'''
Compute maximum identity to other partitions for each sequence
in a partitioning using needleall.
'''
import argparse
import subprocess
import os
import pandas as pd
from tqdm import tqdm
from typing import List, Dict, Tuple
from itertools import groupby
import concurrent.futures
import math

NORMALIZATIONS = {'shortest': lambda a,b,c: a/min(b,c), # a identity b len(seq1) c len(seq2)
                  'longest': lambda a,b,c: a/max(b,c),
                  'mean' : lambda a,b,c: a/((b+c)/2),
                  }


def partition_fasta_file(fasta_file: str, partition_file: str, sep='|') -> Tuple[int, Dict[str, int]]:
    '''
    Make a temporary fasta file for each individual partition, given the .csv assignments.
    Also return the number of partitions and the sequence lengths for convenience, so we don't have to 
    read the file again.
    '''
    seq_dict = {} # Acc : Seq
    seq_lens = {} # Acc : len

    # read fasta
    with open(fasta_file, 'r') as f:

        id_seq_groups = (group for group in groupby(f, lambda line: line.startswith(">")))

        for is_id, id_iter in id_seq_groups:
            if is_id: # Only needed to find first id line, always True thereafter
                id = next(id_iter).strip().split(sep)[0]
                seq = "".join(seq.strip() for seq in next(id_seq_groups)[1])
                seq_dict[id.lstrip('>')] =  seq
                seq_lens[id.lstrip('>')] = len(seq)
    
    # read partition table
    #'AC,label-val,cluster\n'
    df = pd.read_csv(partition_file)
    
    for idx, cl in enumerate(df['cluster'].unique()):
        acc_ids = df.loc[df['cluster'] == cl]['AC']
        with open(f'partition_{idx}.fasta.tmp', 'w') as f:
            for acc_id in acc_ids:
                f.write(f'>{acc_id}\n')
                f.write(seq_dict[acc_id]+'\n')
    return idx + 1, seq_lens


def chunk_fasta_file(fasta_file: str, n_chunks: int, prefix: str = 'chunk', sep: str = '|') -> int:
    '''
    Break up fasta file into multiple smaller files that can be
    used for multiprocessing.
    Returns the number of generated chunks.
    '''

    

    # read fasta
    ids = []
    seqs = []
    with open(fasta_file, 'r') as f:

        id_seq_groups = (group for group in groupby(f, lambda line: line.startswith(">")))
        for is_id, id_iter in id_seq_groups:
            if is_id: # Only needed to find first id line, always True thereafter
                ids.append(next(id_iter).strip().split(sep)[0])
                seqs.append("".join(seq.strip() for seq in next(id_seq_groups)[1]))


    chunk_size = math.ceil(len(ids)/n_chunks)
    empty_chunks = 0
    for i in range(n_chunks):
        # because of ceil() we sometimes make less partitions than specified.
        if i*chunk_size>=len(ids):
            empty_chunks +=1
            continue

        chunk_ids = ids[i*chunk_size:(i+1)*chunk_size]
        chunk_seqs = seqs[i*chunk_size:(i+1)*chunk_size]

        with open(f'{prefix}_{i}.fasta.tmp', 'w') as f:
            for id, seq in zip(chunk_ids, chunk_seqs):
                f.write(id+'\n')
                f.write(seq+'\n')

    return n_chunks - empty_chunks

def compute_edges(query_fp: str,
                  library_fp: str,
                  results_dict: Dict[str, Tuple[float, str]], 
                  seq_lens: Dict[str,int],
                  pbar: tqdm,
                  denominator = 'full',
                  delimiter: str = '|',
                  is_nucleotide: bool = False,
                  gapopen: float = 10,
                  gapextend: float = 0.5,
                  endweight: bool = False,
                  endopen: float = 10,
                  endextend: float = 0.5,
                  matrix: str = 'EBLOSUM62',
                  ) -> None:
    '''
    Run needleall on query_fp and library_fp,
    Retrieve pairwise similiarities, transform and
    insert into edge_dict.
    '''
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
               type_1, type_2, query_fp, library_fp]
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
                
                count += 1
                if count == 1000:
                    pbar.update(1000)
                    count = 0
                # Gaps:           0/142 ( 0.0%)
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
                #line = "# Identity:      14/443 ( 3.2%)"
                # n_matches =  int(line[11:].split('/')[0]) #int() does not mind leading spaces

                if this_qry not in results_dict:
                    results_dict[this_qry] = (identity, this_lib)
                else:
                    if results_dict[this_qry][0]<identity:
                        results_dict[this_qry] = (identity, this_lib)
                








def align_partitions(fasta_file: str, 
                     partition_file: str,
                     denominator = 'full',
                     delimiter: str = '|',
                     is_nucleotide: bool = False,
                     gapopen: float = 10,
                     gapextend: float = 0.5,
                     endweight: bool = False,
                     endopen: float = 10,
                     endextend: float = 0.5,
                     matrix: str = 'EBLOSUM62',
                     n_procs: int = 8) -> Dict[str, Tuple[float, str]]:
    '''
    Align each partition against each other partition,
    and find the shortest distance of each sample to any other sample in 
    any other partition.
    '''
    # Make partition files.
    print('Writing temporary fasta files...')
    n_partitions, seq_lens = partition_fasta_file(fasta_file, partition_file)

    results_dict = {} # Acc: [max_id, Acc]

    part_size = len(seq_lens) // n_partitions
    print('Aligning all partitions...')
    # this is still wrong, we would need to count the actual number of seqs in each partition
    # like that its just an upper bound.
    pbar = tqdm(total=part_size*part_size*(n_partitions*n_partitions - n_partitions)) #inner complexity x nested for loops.
    jobs = []
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_procs) as executor:
        for query_partition in range(n_partitions):
            for lib_partition in range(n_partitions):
                if query_partition == lib_partition:
                    continue
                else:
                    q = f'partition_{query_partition}.fasta.tmp'
                    l = f'partition_{lib_partition}.fasta.tmp'

                    # chunk one of the files to max out threads.
                    # otherwise, can at max use (n_partitions x n_partitions) - n_partitions threads.
                    n_chunks = chunk_fasta_file(q, n_chunks=10, prefix=f'chunk_p_{query_partition}')
                    for i in range(n_chunks):
                        q_c = f'chunk_p_{query_partition}_{i}.fasta.tmp'
                        future = executor.submit(compute_edges, q_c, l, results_dict, seq_lens, pbar, denominator, delimiter, is_nucleotide, gapopen, gapextend, endweight, endopen, endextend, matrix)
                        jobs.append(future)


    for x in os.listdir():
        if 'chunk_p_' in x:
                os.remove(x)

    for i in range(n_partitions):
        os.remove(f'partition_{i}.fasta.tmp')
    return results_dict



def main() -> None:

    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta-file', type=str, help = 'Path to fasta file. needs Graph-Part header format.')
    parser.add_argument('--partition-file', type=str, help='.csv partition assignments')
    parser.add_argument('--out-file', type=str, default = 'homology_reduced_split_with_metrics.csv')

    parser.add_argument('--denominator', type=str, default='full')
    parser.add_argument('--nucleotide', action='store_true')
    parser.add_argument('--gapopen', type=int, default=10)
    parser.add_argument('--gapextend', type=float, default=0.5)
    parser.add_argument('--endweight', action='store_true')
    parser.add_argument('--endopen', type=int, default=10)
    parser.add_argument('--endextend', type=float, default=0.5)
    parser.add_argument('--matrix', type=str, default='EBLOSUM62')
    parser.add_argument('--threads', type=int, default=8)
    args = parser.parse_args()

    max_identities = align_partitions(args.fasta_file, 
                                      args.partition_file, 
                                      args.denominator, 
                                      '|', 
                                      args.nucleotide,
                                      args.gapopen,
                                      args.gapextend,
                                      args.endweight,
                                      args.endopen,
                                      args.endextend,
                                      args.matrix,
                                      args.threads)


    df_partitions = pd.read_csv(args.partition_file)
    df_partitions['max_ident_other'] = None
    df_partitions['closest_other'] = None

    for idx, row in df_partitions.iterrows():

        max_ident, closest = max_identities[row['AC']]
        df_partitions.loc[idx, 'max_ident_other'] = max_ident
        df_partitions.loc[idx, 'closest_other'] = closest

    
    df_partitions.to_csv(args.out_file)

if __name__ == '__main__':
    main()