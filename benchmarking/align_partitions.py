'''
Script to align a fasta file to another. Compute all pairwise identities.
Uses needleall.

Most functions are adapted from needle_utils.py
Result is saved as a hdf5 file with three arrays (seq_a,) (seq_b,) (identity,)
'''

# TODO this whole script sucks - just make a parallel edge list script that dumps to pickle
from typing import Tuple, List, Dict
import argparse
import subprocess
import time
import os
import concurrent.futures
from fasta_utils import parse_fasta, chunk_fasta_file
import h5py


NORMALIZATIONS = {'shortest': lambda a,b,c: a/min(b,c), # a identity b len(seq1) c len(seq2)
                  'longest': lambda a,b,c: a/max(b,c),
                  'mean' : lambda a,b,c: a/((b+c)/2),
                  }

def get_len_dict(ids: List[str], seqs: List[str]) -> Dict[str,int]:
    '''Get a dictionary that contains the length of each sequence.'''
    len_dict = {}
    for id, seq in zip(ids, seqs):
        len_dict[id.lstrip('>')] = len(seq)
    
    return len_dict


def run_needle(query_fp: str,
                library_fp: str,
                identity_dict: Dict[Tuple[str,str],float],
                seq_lens: Dict[str,int],
                denominator = 'full',
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
    command = ["needleall","-auto","-stdout", 
               "-aformat", "pair", 
               "-gapopen", str(gapopen),
               "-gapextend", str(gapextend),
               "-endopen", str(endopen),
               "-endextend", str(endextend),
               "-datafile", matrix,
               "-sprotein1", "-sprotein2", query_fp, library_fp]
    if endweight:
        command = command + ["-endweight"]

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

                
                identity_dict[(this_qry, this_lib)] = identity
                


def run_needle_multithread(fasta_file: str, 
                           ref_file: str, 
                           n_chunks: int, 
                           n_procs: int, 
                           **kwargs) -> Dict[Tuple[str,str],float]:
    '''
    Chunk the input file and start multithreaded alignment.
    Returns identities in a dict {(seqA, seqB): identity}
    '''
    # chunk the input
    ids, seqs = parse_fasta(fasta_file)
    seq_lens = get_len_dict(ids, seqs)
    n_chunks = chunk_fasta_file(ids, seqs, n_chunks) #get the actual number of generated chunks.
    identity_dict = {}
    # start n_procs threads, each thread starts a subprocess
    # Because of threading's GIL we can write edges directly to the full_graph object.
    jobs = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_procs) as executor:
        for i in range(n_chunks):
                q = f'needle_{i}.fasta.tmp'
                future = executor.submit(run_needle, q, ref_file, identity_dict, seq_lens, **kwargs)
                jobs.append(future)

        for job in jobs:
            if job.exception() is not None:
                print(job.exception())
                raise RuntimeError('One of the alignment threads did not complete sucessfully.')


    #delete the chunks
    for i in range(n_chunks):
        os.remove(f'needle_{i}.fasta.tmp')

    return identity_dict





def measure_overlap():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset_name', '-f', type=str)
    parser.add_argument('--dataset_name_2', '-f2', type=str)

    parser.add_argument("-nt","--threads",type=int, help='Number of threads to run in parallel.', default=None)
    parser.add_argument("-nc","--chunks",type=int, help='Number of chunks to split the fasta file.', default=100)

    # customize needle
    parser.add_argument("-dn","--denominator",type=str, help='Denominator to use for sequence identity computation.', 
                        choices=['full', 'shortest', 'longest', 'mean', 'no_gaps'], 
                        default='full',
                        )
    parser.add_argument('--gapopen','-gapopen', type=float, default=10, help='Passed to needle. See EMBOSS documentation.')
    parser.add_argument('--gapextend','-gapextend', type=float, default=0.5, help='Passed to needle. See EMBOSS documentation.')
    parser.add_argument('--endweight','-endweight', action='store_true', help='Passed to needle. See EMBOSS documentation.')
    parser.add_argument('--endopen','-endopen', type=float, default=10, help='Passed to needle. See EMBOSS documentation.')
    parser.add_argument('--endextend','-endextend', type=float, default=10, help='Passed to needle. See EMBOSS documentation.')
    parser.add_argument('--matrix', '--datafile','-datafile', type=str, default='EBLOSUM62', help='Passed to needle. See EMBOSS documentation.')


    args = parser.parse_args()

    # 2. get partition files
    name_1 = os.path.basename(args.dataset_name)
    name_2 = os.path.basename(args.dataset_name_2)

    print(f'Using {args.threads} threads to align {name_1} to {name_2}')
    print(f'Processing {name_1} as {args.chunks} chunks.')
    identity_dict = run_needle_multithread(
                                        args.dataset_name, 
                                        args.dataset_name_2, 
                                        n_chunks=args.chunks, 
                                        n_procs=args.threads,
                                        denominator = args.denominator,
                                        gapopen = args.gapopen,
                                        gapextend=args.gapextend,
                                        endweight=args.endweight,
                                        endopen=args.endopen,
                                        endextend=args.endextend,
                                        matrix=args.matrix
                                        )

    dump_name = f'{name_1[:-6]}_to_{name_2[:-6]}.hdf5'

    f = h5py.File(dump_name, 'w')

    dset_idents = f.create_dataset("identities", (len(identity_dict),), dtype=float)
    dset_qry = f.create_dataset("query", (len(identity_dict),), dtype="S10")
    dset_lib = f.create_dataset("library", (len(identity_dict),), dtype="S10")

    for idx, (key, value) in enumerate(identity_dict.items()):
        qry, lib = key
        dset_idents[idx] = value
        dset_qry[idx] = qry.encode("ascii", "ignore")
        dset_lib[idx] = lib.encode("ascii", "ignore")

    f.close()  

# asciiList = [n.encode("ascii", "ignore") for n in strList]
# h5File.create_dataset('xxx', (len(asciiList),1),'S10', asciiList)
    #save_path = name[:-6] + '_num_entities_above.json'
    #tabulate_percent_overlap(id_dict, save_path)

    # visualizations: grouped boxplots of cross-identities per partition, red line for threshold
    # table with % above 
    # % above defined as % proteins with at least 1 identity above threshold.




if __name__ == '__main__':
    s = time.perf_counter()
    measure_overlap()
    elapsed = time.perf_counter() - s
    print(f"{__file__} executed in {elapsed:0.2f} seconds.")