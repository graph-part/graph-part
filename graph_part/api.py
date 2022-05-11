'''
Python interface for Graph-Part.
'''
from typing import Iterable, List, Dict, Union, Tuple
import numpy as np
import os
import pandas as pd
from .graph_part import run_partitioning


def _write_fasta(sequences: Dict[str,str], labels: Dict[str,str] = None, priority: Dict[str,str] = None, fname='graphpart_api.fasta.tmp') -> None:
    '''Make a fasta file from the provided dictionaries. This will be used for alignment.'''

    # No IDs provided-make our own.
    with open(fname, 'w') as f:
        for name, seq in sequences.items():
            header = f'>{name}'
            if labels is not None:
                header = header + f'|label={labels[name]}'
            if priority is not None:
                header = header + f'|priority={priority[name]}'

            f.write(header + '\n')
            f.write(seq + '\n')
            

def _convert_to_dict(sequences: Union[List[str], np.ndarray, pd.core.series.Series, Dict[str,str]],
                    labels: Union[List[str], np.ndarray, pd.core.series.Series, Dict[str,str]] = None,
                    priority: Union[List[str], np.ndarray, pd.core.series.Series, Dict[str,str]] = None,
                    ) -> Tuple[Dict[str,str], Union[None, Dict[str,str]], Union[None, Dict[str,str]]]:
    '''
    For simplicity, we process all input data as dicts internally. Do not allow user to mix
    input types between dicts and arrays/series.
    '''

    # ensure that there is no dict/list mixup in the inputs. Could handle interally,
    # but probably better to force the user to avoid spurious mismatch errors.
    if type(sequences) == dict:
        if labels is not None and type(labels) != dict:
            raise ValueError('When providing sequences as dictionary, also provide labels as dictionary.')
        if priority is not None and type(priority) != dict:
            raise ValueError('When providing sequences as dictionary, also provide priorities as dictionary.')

        return sequences, labels, priority
    
    if type(sequences) == pd.core.series.Series:
        if labels is not None and type(labels) !=  pd.core.series.Series:
            raise ValueError('When providing sequences as series, also provide labels as series.')
        if priority is not None and type(priority) !=  pd.core.series.Series:
            raise ValueError('When providing sequences as series, also provide priorities as series.')

        return sequences.to_dict(), labels.to_dict(), priority.to_dict()

    elif type(sequences) in [np.ndarray, pd.core.series.Series, list]:
        if labels is not None and type(labels) not in [np.ndarray,  list]:
            raise ValueError('When providing sequences as array, also provide labels as array.')
        if priority is not None and  type(priority) not in [np.ndarray, list]:
            raise ValueError('When providing sequences as array, also provide priorities as array.')

        seq_dict = {}
        lab_dict = {} if labels is not None else None
        pri_dict = {} if priority is not None else None
        for i in range(len(sequences)):
            seq_dict[f'seq_{i}'] = sequences[i]
            if labels is not None:
                lab_dict[f'seq_{i}'] = labels[i]
            if priority is not None:
                pri_dict[f'seq_{i}'] = priority[i]

        return seq_dict, lab_dict, pri_dict



def stratified_k_fold(sequences: Union[List[str], np.ndarray, Dict[str,str]], 
                     labels: Union[List[str], np.ndarray, Dict[str,str]] = None,
                     priority: Union[List[str], np.ndarray, Dict[str,str]] = None,
                     partitions: int = 5,
                     threshold: float = 0.3,
                     transformation: str = 'one-minus',
                     alignment_mode: str = 'mmseqs2',
                     initialization_mode: str = 'slow-nn',
                     no_moving: bool = False,
                     remove_same: bool = False,
                     save_checkpoint_path: str = None,
                     denominator: str = 'full',
                     nucleotide: bool = False,
                     triangular: bool = False,
                     threads: int = 4,
                     chunks: int = 10,
                     parallel_mode: str = 'multithread',
                     gapopen: float = 10,
                     gapextend: float = 0.5,
                     endweight: bool = False,
                     endopen: float =10,
                     endextend: float = 0.5,
                     matrix: str = 'EBLOSUM62',
                     edge_file: str = None,
                     metric_column: str = None,
                     ) -> List[Iterable]:
    '''
    Split an array or dictionary of sequences into balanced k folds.

    Parameters:
    --------
        sequences : list, dict, numpy array or pandas series.
            Iterable of protein or nucleotide sequences to partition.

        labels :  list, dict, numpy array or pandas series.
            Optional iterable of labels.

        priority : list, dict, numpy array or pandas series.
            Optional iterable of priority labels.

        threshold : float
            Percent identity threshold for paritioning.
    

    Returns:
    ---------
        splitting: list, length = partitions
            List of ids belonging to each partition. If the input was an array or list, this will contain indices, else the sequence IDs.

    #TODO add warnings that arguments will be ignored depending on alignment_mode.
    '''
    # 1. Validate arguments.
    if alignment_mode not in ['mmseqs2', 'needle', 'precomputed']:
        raise NotImplementedError(f'Alignment mode {alignment_mode} is not implemented. Choose either `needle` or `mmseqs2`.')

    # Sort out the input formats.
    original_type = type(sequences)
    sequences, labels, priority = _convert_to_dict(sequences, labels, priority)

    # 2. Write the data to a temporary fasta file for alignment.
    _write_fasta(sequences, labels, priority, 'graphpart_api.fasta.tmp')

    config = {
        "alignment_mode": alignment_mode,
        "fasta_file": "graphpart_api.fasta.tmp",
        "threshold": threshold,
        "partitions": partitions,
        "transformation": transformation,
        "out_file": "graphpart_python", # not used anyway
        "priority_name": "priority" if priority is not None else None,
        "labels_name": "label" if labels is not None else None,
        "initialization_mode": initialization_mode,
        "no_moving": no_moving,
        "remove_same": remove_same,
        "test_ratio": 0,
        "val_ratio": 0,
        "save_checkpoint_path": save_checkpoint_path,
        "denominator": denominator,
        "nucleotide": nucleotide,
        "triangular": triangular,
        "threads": threads,
        "chunks": chunks,
        "parallel_mode": parallel_mode,
        "gapopen": gapopen,
        "gapextend": gapextend,
        "endweight": endweight,
        "endopen": endopen,
        "endextend": endextend,
        "matrix": matrix,
        "edge_file": edge_file,
        "metric_column": metric_column,
        "allow_moving": not no_moving, # silly conversions because in the CLI we want to have those default-false.
        "removal_type": not remove_same,
    }

    # 3. Partition
    partition_assignment_df = run_partitioning(config, write_output_file=False, write_json_report=False, verbose=False)
    os.remove(config['fasta_file'])

    # 4. Make output lists.
    partition_assignment_df = partition_assignment_df.reset_index()
    outs = []
    # iterate over all created folds and add to the output list.
    for _, sub_df in partition_assignment_df.groupby('cluster'):

        if original_type in [np.ndarray, list]:
            outs.append(sub_df.index.tolist())
        else:
            outs.append(sub_df['AC'].tolist())
    return outs




def train_test_validation_split(sequences: Union[List[str], np.ndarray, Dict[str,str]], 
                     labels: Union[List[str], np.ndarray, Dict[str,str]] = None,
                     priority: Union[List[str], np.ndarray, Dict[str,str]] = None,
                     test_size: float = 0.15, 
                     valid_size: float = 0,
                     threshold: float = 0.3,
                     transformation: str = 'one-minus',
                     alignment_mode: str = 'mmseqs2',
                     initialization_mode: str = 'slow-nn',
                     no_moving: bool = False,
                     remove_same: bool = False,
                     save_checkpoint_path: str = None,
                     denominator: str = 'full',
                     nucleotide: bool = False,
                     prefilter: bool = False,
                     triangular: bool = False,
                     threads: int = 4,
                     chunks: int = 10,
                     parallel_mode: str = 'multithread',
                     gapopen: float = 10,
                     gapextend: float = 0.5,
                     endweight: bool = False,
                     endopen: float =10,
                     endextend: float = 0.5,
                     matrix: str = 'EBLOSUM62',
                     edge_file: str = None,
                     metric_column: str = None,
                     ) -> List[Iterable]:
    '''
    Split an array or dictionary of sequences into train-validation-test subsets.

    Parameters:
    --------
        sequences : list, dict, numpy array or pandas series.
            Iterable of protein or nucleotide sequences to partition.

        labels :  list, dict, numpy array or pandas series.
            Optional iterable of labels.

        priority : list, dict, numpy array or pandas series.
            Optional iterable of priority labels.

        threshold : float
            Percent identity threshold for paritioning.
    

    Returns:
    ---------
        splitting: list, length = 2 or 3 if valid_size>0
            List of ids belonging to each partition. If the input was an array or list, this will contain indices, else the sequence IDs.

    #TODO add warnings that arguments will be ignored depending on alignment_mode.

    '''
    # 1. Validate arguments.
    if alignment_mode not in ['mmseqs2', 'needle', 'precomputed']:
        raise NotImplementedError(f'Alignment mode {alignment_mode} is not implemented. Choose either `needle` or `mmseqs2`.')

    if test_size*100 % 5 !=0 or valid_size*100 % 5 != 0:
        raise NotImplementedError('Graph-Part currently only supports ratios that are a multiple of 0.05!')

    if test_size* 100 % 10 == 0 and valid_size* 100 % 10 == 0:
        partitions = 10
    else:
        partitions = 20

    # if test_ratio is 0 but val_ratio is defined, just swap the two. Then everything works.
    if test_size ==0:
        test_size = valid_size
        valid_size = test_size



    # Sort out the input formats.
    original_type = type(sequences)
    sequences, labels, priority = _convert_to_dict(sequences, labels, priority)

    # 2. Write the data to a temporary fasta file for alignment.
    _write_fasta(sequences, labels, priority, 'graphpart_api.fasta.tmp')

    config = {
        "alignment_mode": alignment_mode,
        "fasta_file": "graphpart_api.fasta.tmp",
        "threshold": threshold,
        "partitions": partitions,
        "transformation": transformation,
        "out_file": "graphpart_python", # not used anyway
        "priority_name": "priority" if priority is not None else None,
        "labels_name": "label" if labels is not None else None,
        "initialization_mode": initialization_mode,
        "no_moving": no_moving,
        "remove_same": remove_same,
        "test_ratio": test_size,
        "val_ratio": valid_size,
        "save_checkpoint_path": save_checkpoint_path,
        "denominator": denominator,
        "nucleotide": nucleotide,
        "prefilter": prefilter,
        "triangular": triangular,
        "threads": threads,
        "chunks": chunks,
        "parallel_mode": parallel_mode,
        "gapopen": gapopen,
        "gapextend": gapextend,
        "endweight": endweight,
        "endopen": endopen,
        "endextend": endextend,
        "matrix": matrix,
        "edge_file": edge_file,
        "metric_column": metric_column,
        "allow_moving": not no_moving, # silly conversions because in the CLI we want to have those default-false.
        "removal_type": not remove_same,
    }

    # 3. Partition
    partition_assignment_df = run_partitioning(config, write_output_file=False, write_json_report=False, verbose=False)
    os.remove(config['fasta_file'])

    # 4. Make output lists.
    partition_assignment_df = partition_assignment_df.reset_index()
    outs = []
    # iterate over all created folds and add to the output list.
    for _, sub_df in partition_assignment_df.groupby('cluster'):

        if original_type in [np.ndarray, list]:
            outs.append(sub_df.index.tolist())
        else:
            outs.append(sub_df['AC'].tolist())
    return outs

