import argparse
import subprocess
import os
import re
import shutil
from typing import List, Tuple, Dict
import pathlib
import numpy as np



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
    for i in u_cluster:
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


# reduction_result.clstr
def parse_clustering(fpath) -> Tuple[List[str], List[str]]:
    '''
    Parse CD-HIT clustering output.
    Does not return full headers, so we only retrieve the ACC of each seq.
    '''
    clusters = []
    accs = []

    current_cluster = -1
    with open(fpath) as f:
        for line in f:
            if line.startswith('>'):
                current_cluster += 1
            else:
                spl = line.strip().split(' ')
                for s in spl:
                    if s.startswith('>'):
                        acc = s.split('|')[0].lstrip('>')

                clusters.append(current_cluster)
                accs.append(acc)
            
    
    return clusters, accs


def get_labels(fasta_file: str, labels_name: str = 'label') -> Dict[str, str]:
    
    labels = {}
    label_id_dict = {}
    label_count = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                spl = line.strip().split('|')
                acc = spl[0].lstrip('>')
                label = '0'
                for s in spl[1:]:
                    if '=' in s:
                        param_spl = s.split('=')
                        if param_spl[0] == labels_name:
                            label = str(param_spl[1].strip())

                # convert to integers.
                if label not in label_id_dict:
                    label_id_dict[label]  = label_count
                    label_count += 1
                
                labels[acc] = label_id_dict[label]
            

    return labels


def psicdhit_homology_cluster(entity_fp: str, threshold: float = 0.3) -> Tuple[List[int], List[str]]:

    out_prefix = 'reduction_result'
   
    subprocess.run(['cd-hit', '-i', entity_fp, '-o', 'pre_reduced_1' , '-c', '0.9', '-n', '5', '-T', '0' ])
    subprocess.run(['cd-hit', '-i', 'pre_reduced_1', '-o', 'pre_reduced_2' , '-c', '0.6', '-n', '4', '-T', '0' ]) #1792

    psi_path = pathlib.Path(__file__).parent.resolve() / 'psi_cd_hit.pl'
        #subprocess.run(['perl','-I', psi_path.parent, psi_path,
    subprocess.run(['perl', psi_path, '-i', 'pre_reduced_2', '-o', out_prefix, '-c', str(threshold)])


    clusters, accs = parse_clustering('reduction_result.clstr')
    
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

    return clusters, accs



def cdhit_homology_cluster(entity_fp: str, threshold: float = 0.3, is_nucleotide: bool = False) -> Tuple[List[int], List[str]]:

    if not is_nucleotide and threshold<0.4:
        return psicdhit_homology_cluster(entity_fp, threshold)
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

    clusters, accs = parse_clustering('reduction_result.clstr')
    
    for x in os.listdir():
        if 'reduction_result' in x:
            if os.path.isdir(x):
                shutil.rmtree(x)
            else:
                os.remove(x)
    
    return clusters, accs



def main() -> None:

    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta-file', type=str, help = 'Path to fasta file. needs Graph-Part header format.')
    parser.add_argument('--out-file', type=str, default = 'homology_reduced_split.csv')
    parser.add_argument('--labels-name', default=None)
    parser.add_argument('-th', '--threshold', type=float, default = 0.3)
    parser.add_argument('-pa', '--partitions', type=int, default=5)

    args = parser.parse_args()

    cluster_ids, acc_ids =  cdhit_homology_cluster(args.fasta_file, args.threshold)

    labels_dict =  get_labels(args.fasta_file, args.labels_name)

    labels = [labels_dict[x] for x in acc_ids]
    

    # partition the data.
    cluster_ids = np.array(cluster_ids)
    labels = np.array(labels)
    dummy_kingdoms = np.zeros_like(labels)
    cl_number = partition_assignment(cluster_ids, dummy_kingdoms, labels, args.partitions, len(np.unique(labels)), 1)



    with open(args.out_file, 'w') as f:
        f.write('AC,label-val,cluster\n')

        for idx, acc in enumerate(acc_ids):
            f.write(f'{acc},{labels[idx]},{cl_number[idx]}\n')





if __name__ == '__main__':
    main()
