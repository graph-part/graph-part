'''
Code to generate train-val-test splits using Graph-Part.
We first partition the data into 10 or 20 partitions. These
are then combined to yield the splits. Then removal is applied,
as always.
'''
import networkx as nx
import numpy as np
import pandas as pd
from collections import Counter
from itertools import combinations
from typing import List, Tuple


def check_train_val_test_args(args):
    '''If a train-val-test split is to be done, we
    need to fix the "partitions" to be compatible.'''
    if args.test_ratio*100 % 5 !=0 or args.val_ratio*100 % 5 != 0:
        raise NotImplementedError('Graph-Part currently only supports ratios that are a multiple of 0.05!')

    if args.test_ratio* 100 % 10 == 0 and args.val_ratio* 100 % 10 == 0:
        setattr(args, 'partitions', 10)
    else:
        setattr(args, 'partitions', 20)

    # if test_ratio is 0 but val_ratio is defined, just swap the two. Then everything works.
    if args.test_ratio ==0:
        setattr(args, 'test_ratio', args.val_ratio)
        setattr(args, 'val_ratio', 0.0)
    
    
    

def compute_partition_similarity_matrix(full_graph: nx.classes.graph.Graph, part_graph: nx.classes.graph.Graph, n_partitions: int, threshold: float) -> np.ndarray:
    '''Compute a similarity matrix of the partitions. Metric = number of connections between.'''
    partition_connections = np.zeros((n_partitions, n_partitions))
    #iterate over all sequences
    for n,d in full_graph.nodes(data=True):
        # get the partition of the sequence 
        self_cluster = part_graph.nodes[n]['cluster']
        #get the neighbors of the sequence
        neighbours = nx.neighbors(full_graph,n) 
        #count the number of neighbors of each partition
        neighbour_clusters = Counter((part_graph.nodes[nb]['cluster'] for nb in nx.neighbors(full_graph,n) if full_graph[n][nb]['metric'] < threshold))

        for cl, count in neighbour_clusters.items():
            partition_connections[int(self_cluster), int(cl)] += count #Partition id ['cluster'] is float.

    return partition_connections

def find_best_partition_combinations(partition_connections: np.ndarray, n_train: int, n_test: int) -> Tuple[List[int], List[int], List[int]]:
    '''
    Brute force try all combinations of partitions to find the set of partitions that have the maximum connections to each other.
    This works because the expected number of partitions is low enough, e.g. steps of 10% or 5% -> max 20 partitions.
    Don't do this when using 1% steps, will explode.
    '''
    partitions = list(range(partition_connections.shape[0]))
    
    
    def get_best_combination(partitions, n):
        best_combination = []
        best_score = -1 # When a partition has no connection to any other, or to itself, the score will be 0. So make this -1
        for comb in combinations(partitions, n):
            score = partition_connections[comb,:][:,comb].sum()
            if score>best_score:
                best_combination=comb
                best_score = score
                
        return best_combination
    
    # find the best combination for train.
    train_partitions = get_best_combination(partitions, n_train)
    remainder = [x for x in partitions if x not in train_partitions]
    # find the best combination for test.
    test_partitions = get_best_combination(remainder, n_test)
    
    # the remainder is the validation set.
    val_partitions = [x for x in remainder if x not in test_partitions]
    
    return train_partitions, test_partitions, val_partitions


def train_val_test_split(part_graph: nx.classes.graph.Graph, 
                     full_graph: nx.classes.graph.Graph, 
                     threshold: float, 
                     test_ratio: float,
                     val_ratio: float, 
                     n_partitions: int = 10) -> None:
    '''
    Merge pre-removal partitions to generate a train-val-test split.
    '''
    n_train = int(round(n_partitions * (1-val_ratio-test_ratio))) # has a .999999999 float issue without rounding.
    n_test = int(n_partitions * test_ratio)
    n_val = int(n_partitions * val_ratio)

    # For each partition, measure the overlap to other partitions.
    # partition_connections is essentially a similarity matrix of all the partitions.
    partition_connections = compute_partition_similarity_matrix(full_graph, part_graph, n_partitions, threshold)

    # Given the similarity matrix, find the combinations with maximum overlap.
    # By doing this now, we reduce the number of move/removal operations later.
    train_partitions, test_partitions, val_partitions = find_best_partition_combinations(partition_connections, n_train, n_test)
    # Given the new assignments, update the graph.
    # part_graph has the following format: {'C0IW58': {'cluster': 0.0, 'C-size': 124, 'label-counts': array([92, 32])}

    # Compute the new statistics to add to the nodes.
    df = pd.DataFrame(((d) for n,d in full_graph.nodes(data=True)))
    df['cluster'] = [part_graph.nodes[n]['cluster'] for n in full_graph.nodes()]
    df['AC'] = [n for n in full_graph.nodes()]
    df = df.groupby(['cluster','label-val'])['AC'].count().reset_index().pivot_table(values='AC',columns=['label-val'],index=['cluster'])

    train_statistics = df.loc[[float(x) for x in train_partitions]].sum(axis=0)
    train_attributes = {'cluster': 0.0, 'C-size': sum(train_statistics), 'label-counts': train_statistics.to_numpy()}
    test_statistics = df.loc[[float(x) for x in test_partitions]].sum(axis=0)
    test_attributes = {'cluster': 1.0, 'C-size': sum(test_statistics), 'label-counts': test_statistics.to_numpy()}
    val_statistics = df.loc[[float(x) for x in val_partitions]].sum(axis=0)
    val_attributes = {'cluster': 2.0, 'C-size': sum(val_statistics), 'label-counts': val_statistics.to_numpy()}


    # Now, update part_graph
    for n, d in part_graph.nodes(data=True):
        if int(d['cluster']) in train_partitions:
            nx.set_node_attributes(part_graph, {n:train_attributes})
        elif int(d['cluster']) in test_partitions:
            nx.set_node_attributes(part_graph, {n:test_attributes})
        else:
            nx.set_node_attributes(part_graph, {n:val_attributes})



# Not used.
def find_best_partition_combinations_heuristic(partition_connections: np.ndarray, n_train: int, n_test: int):
    '''
    Brute force try combinations of partitions to find the set of partitions that have the maximum connections to each other.
    This here still tries a lot, but is not exhaustive. 
    Version above: O(min(n^k, n^(n-k)))
    This:        : O(n*n*k)
    As long as we have 10 or 20 partitions, n choose k is fine. This works too, but untested for larger n.
    '''
    # 1. Get train partitions
    best_combination = []
    best_score = 0
    current_score = 0
    
    # Try each partition as the starting point. complexity n*n*k
    for partition_id in range(partition_connections.shape[0]): # complexity n
        current_combination = []
        current_combination.append(partition_id)
        
        # Add partitions to the set until the size is reached.
        while(len(current_combination)<n_train): #complexity k
            
            # find the best partition to add.
            best_connections = 0
            best_connected = None
            for i in range(partition_connections.shape[0]): # complexity n
                #print('Current', current_combination, 'testing', i, 'best', best_connections)
                if i in current_combination:
                   # print('i in current')
                    continue
                
                connections = partition_connections[i, current_combination].sum()
                #print(f'Connections for {i}: {connections}')
                if connections > best_connections:
                    best_connected = i
                    best_connections = connections
            
            current_combination.append(best_connected)
        
        # Got a set. Keep if best, else discard.
        # We double count here (do np.triu to avoid), but it does not matter if we do it for all.
        current_score = partition_connections[current_combination,:][:,current_combination].sum()
        if current_score>best_score:
            best_combination = current_combination
            best_score = current_score
        
    return best_combination