import argparse
import subprocess
import os
import re
import shutil
from typing import List, Tuple, Dict
import pathlib
import numpy as np
from tqdm.auto import tqdm
import time
import json

def partition_assignment(cluster_vector : np.array, kingdom_vector: np.array, label_vector: np.array, n_partitions: int, n_class: int, n_kingdoms: int) -> np.array:
    ''' Function to separate proteins into N partitions with balanced classes
        Inputs:
            cluster_vector: (n_sequences) integer vector of the cluster id of each sequence
            kingdom_vector: (n_sequences) integer vector of the taxonomy group of each sequence
            label_vector  : (n_sequences) integer vector of another balancing criterion
        Returns:
            (n_sequences) split indicator vector
    '''
    
    # Unique cluster number
    u_cluster = np.unique(cluster_vector)
    
    # Initialize matrices
    loc_number = np.ones((n_partitions,n_kingdoms,n_class))
    cl_number = np.zeros(cluster_vector.shape[0])
    
    processed_clusters = 1
    n_clusters = u_cluster.shape[0]
    for i in tqdm(u_cluster, desc='Assigning clusters'):
        #if i % 100 == 0:
        #    logger.info(f'Processing cluster {processed_clusters}/{n_clusters}')
        processed_clusters +=1
        # Extract the labels for the proteins in that cluster
        positions = np.where(cluster_vector == i)
        cl_labels = label_vector[positions]
        cl_kingdom = kingdom_vector[positions]
        cl_all = np.column_stack((cl_kingdom,cl_labels))
 
        # Count number of each class
        u, count = np.unique(cl_all, axis=0, return_counts=True)
        
        temp_loc_number = np.copy(loc_number)
        temp_loc_number[:,u[:,0],u[:,1]] += count
        loc_per = loc_number/temp_loc_number
        best_group = np.argmin(np.sum(np.sum(loc_per, axis=2),axis=1))
        loc_number[best_group,u[:,0],u[:,1]] += count
        
        # Store the selected partition
        cl_number[positions] = best_group
    
    print(loc_number.astype(np.int32)-np.ones(loc_number.shape))
    return cl_number


def mmseqs2_homology_cluster(entity_fp: str, threshold: float = 0.3, cluster_mode: int = 0, seq_id_mode: int = 0) -> Tuple[List[int], List[str]]:

    out_prefix = 'reduction_result'
    temp_path = 'mmseqs_temp'
    os.makedirs(temp_path, exist_ok=True)

    subprocess.run(['mmseqs', 'easy-cluster', '--min-seq-id', str(threshold), '--seq-id-mode', str(seq_id_mode), '--cluster-mode', str(cluster_mode), entity_fp, out_prefix, temp_path])
    shutil.rmtree(temp_path)

    accs = []
    clusters = []
    current_cluster = -1
    last_cluster = None
    with open(f'{out_prefix}_cluster.tsv') as f:
        for line_nr, line in enumerate(f):
            spl = line.strip().split('\t')

            rep = spl[0]
            acc = spl[1]

            if rep != last_cluster:
                last_cluster = rep
                current_cluster +=1
            
            accs.append(acc)
            clusters.append(rep)

    for x in os.listdir():
        if 'reduction_result' in x:
            os.remove(x)

    return clusters, accs

def get_labels(identifiers: List[str], labels_name: str = 'label') -> List[str]:

    label_ids = []
    label_dict = {}
    label_count = 0

    for ident in identifiers:

        spl = ident.strip().split('|')
        label = '0'
        for s in spl[1:]:
            if '=' in s:
                param_spl = s.split('=')
                if param_spl[0] == labels_name:
                    label = str(param_spl[1].strip())

        if label not in label_dict:
            label_dict[label]  = label_count
            label_count += 1
        
        label_ids.append(label_dict[label])

    return label_ids


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

    cluster_ids, identifiers =  mmseqs2_homology_cluster(args.fasta_file, args.threshold, args.cluster_mode, args.seq_id_mode)
    out_dict['time_clustering_complete'] = time.perf_counter()
    
    labels =  get_labels(identifiers, args.labels_name)
    acc_ids = [x.split('|')[0] for x in identifiers]

    
    # partition the data.
    cluster_ids = np.array(cluster_ids)
    labels = np.array(labels)
    dummy_kingdoms = np.zeros_like(labels)
    cl_number = partition_assignment(cluster_ids, dummy_kingdoms, labels, args.partitions, len(np.unique(labels)), 1)


    with open(args.out_file, 'w') as f:
        f.write('AC,label-val,cluster\n')

        for idx, acc in enumerate(acc_ids):
            f.write(f'{acc},{labels[idx]},{cl_number[idx]}\n')

    out_dict['time_script_complete'] = time.perf_counter()
    json.dump(out_dict, open(os.path.splitext(args.out_file)[0]+'_report.json','w'))




if __name__ == '__main__':
    main()
