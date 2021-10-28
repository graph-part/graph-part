'''
Convenience script to perform all-vs-all
alignment of two files using needleall.
Produces a .csv file that list the maximum identity in the other file
for each sequence.
'''
import argparse
import subprocess
import os
from tqdm import tqdm
import pandas as pd
import math
from typing import List, Dict, Tuple
from itertools import groupby
import concurrent.futures


NORMALIZATIONS = {'shortest': lambda a,b,c: a/min(b,c), # a identity b len(seq1) c len(seq2)
                  'longest': lambda a,b,c: a/max(b,c),
                  'mean' : lambda a,b,c: a/((b+c)/2),
                  }


def parse_fasta(fastafile: str, sep='|') -> Tuple[List[str],List[str]]:
	'''
    Parses fasta file into lists of identifiers and sequences.
	Can handle multi-line sequences and empty lines.
    Needleall seems to fail when a '|' is between the identifier and the rest of
    the fasta header, so we split the identifier and only return that.
    '''
	ids = []
	seqs = []
	with open(fastafile, 'r') as f:

		id_seq_groups = (group for group in groupby(f, lambda line: line.startswith(">")))

		for is_id, id_iter in id_seq_groups:
			if is_id: # Only needed to find first id line, always True thereafter
			    ids.append(next(id_iter).strip().split(sep)[0])
			    seqs.append("".join(seq.strip() for seq in next(id_seq_groups)[1]))
        
	return ids, seqs


def chunk_fasta_file(ids: List[str], seqs: List[str], n_chunks: int, prefix: str = 'graphpart') -> int:
    '''
    Break up fasta file into multiple smaller files that can be
    used for multiprocessing.
    Returns the number of generated chunks.
    '''

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

def get_len_dict(ids: List[str], seqs: List[str]) -> Dict[str,int]:
    '''Get a dictionary that contains the length of each sequence.'''
    len_dict = {}
    for id, seq in zip(ids, seqs):
        len_dict[id.lstrip('>')] = len(seq)
    
    return len_dict


def compute_edges(query_fp: str,
                  library_fp: str,
                  progress_bar: tqdm,
                  results_dict: Dict[str, Tuple[float, str]], 
                  seq_lens: Dict[str,int],
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
                progress_bar.update(1)

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

                if this_lib not in results_dict:
                    results_dict[this_lib] = (identity, this_qry)
                else:
                    if results_dict[this_lib][0]<identity:
                        results_dict[this_lib] = (identity, this_qry)








def align_partitions(fasta_file_1: str, 
                     fasta_file_2: str,
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
    # Chunk the files.
    ids_1, seqs_1 = parse_fasta(fasta_file_1)
    ids_1 =  [x + '_file_1' for x in ids_1]
    ids_2, seqs_2 = parse_fasta(fasta_file_2)
    ids_2 =  [x + '_file_2' for x in ids_2]

    n_chunks_1 = chunk_fasta_file(ids_1, seqs_1, 100, prefix='file1')
    n_chunks_2 = chunk_fasta_file(ids_2, seqs_2, 100, prefix='file2')
    
    seq_lens = get_len_dict(ids_1, seqs_1)
    seq_lens = seq_lens.update(get_len_dict(ids_2, seqs_2))


    # Align all seqs in file 1 against all seqs in file 2.

    results_dict = {} # Acc: [max_id, Acc]

    jobs = []
    pbar = tqdm(total= len(ids_1) * len(ids_2))
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_procs) as executor:
        for query_chunk in range(n_chunks_1):
            for lib_chunk in range(n_chunks_2):
                
                # chunk query.
                q = f'file1_{query_chunk}.fasta.tmp'
                l = f'file2_{lib_chunk}.fasta.tmp'
                future = executor.submit(compute_edges, q, l, pbar, results_dict, seq_lens, denominator, delimiter, is_nucleotide, gapopen, gapextend, endweight, endopen, endextend, matrix)
                jobs.append(future)



    for i in range(n_chunks_1):
        os.remove(f'file1_{i}.fasta.tmp')
    for i in range(n_chunks_2):
        os.remove(f'file2_{i}.fasta.tmp')
    return results_dict



def main() -> None:

    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta-file-1', type=str, help = 'Path to fasta file. needs Graph-Part header format.')
    parser.add_argument('--fasta-file-2', type=str, help = 'Path to fasta file. needs Graph-Part header format.')
    parser.add_argument('--out-file', type=str, default = 'max_overlaps.csv')

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

    max_identities = align_partitions(args.fasta_file_1, 
                                      args.fasta_file_2, 
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


    
    df_partitions = pd.DataFrame.from_dict(max_identities, columns=['max_ident_other', 'closest_other'], orient='index')
    df_partitions.index.name ='AC'
    #df_partitions = 
    # df_partitions['max_ident_other'] = None
    # df_partitions['closest_other'] = None

    # for idx, row in df_partitions.iterrows():

    #     max_ident, closest = max_identities[row['AC']]
    #     df_partitions.loc[idx, 'max_ident_other'] = max_ident
    #     df_partitions.loc[idx, 'closest_other'] = closest

    
    df_partitions.to_csv(args.out_file)

if __name__ == '__main__':
    main()