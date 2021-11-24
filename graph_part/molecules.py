'''
Functions to run Graph-Part on small molecule data.
We take molecules as SMILE strings, as then we can reuse
a lot of the sequence-based code.
'''
from typing import List, Dict, Tuple, Union, Iterable
import numpy as np
import pandas as pd
import networkx as nx
from tqdm.auto import tqdm
from .graph_part import partition_and_remove


# TODO
# accept pre-converted RDkit molecules, not just smiles.
# additional fingerprint algorithms
# parallelize BulkTanimotoSimilarity
# TODO rdkit tanimoto similarity workflow

def load_entities(molecules: Dict[str,str], labels: Dict[str,str] = None, priorities: Dict[str,str] = None):
    '''
    Construct the graphs by adding nodes. No edges are generated in this step.
    '''
    part_graph = nx.Graph()
    full_graph = nx.Graph()

    labels_out = {}
    for id, mol in molecules.items():
        
        label = '0'
        priority = False
        if labels is not None:
            label = labels[id]
        if priorities is not None:
            priority = priorities[id]

        if label not in labels_out:
            labels_out[label] = {'val':len(labels_out), 'num':0}
        
        labels_out[label]['num'] += 1
        label = labels_out[label]['val']

        full_graph.add_node(id)
        part_graph.add_node(id)
        node_data = {
        'priority': priority,
        'label-val': label
        }
        nx.set_node_attributes(full_graph, {id:node_data})

            
    return full_graph, part_graph, labels_out




def compute_fingerprint_tanimoto_distances(full_graph: nx.classes.graph.Graph, molecules: Dict[str, str], threshold: float) -> None:
    '''Compute the fingerprint tanimoto distances of all pairs and add edges to the graph if they are below
    the threshold'''
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem

    # make the fingerprints
    names = list(molecules.keys())
    mols  = [Chem.MolFromSmiles(x) for x in molecules.values()]

    # TODO make choice of fingerprint algo controllable
    fps = [AllChem.GetMorganFingerprintAsBitVect(x, 2, 1024) for x in mols]

    #ms = [Chem.MolFromSmiles('CCOC'), Chem.MolFromSmiles('CCO'), Chem.MolFromSmiles('COC')]
    #fps = 
    for idx_1 in tqdm(range(len(names))):
        name_1 = names[idx_1]

        similarities = DataStructs.BulkTanimotoSimilarity(fps[idx_1], fps[idx_1+1:])

        for idx_2 in range(len(similarities)):
            name_2 = names[idx_1+1+idx_2]


            metric = 1 - similarities[idx_2]
            if metric > threshold:
                continue
            if not full_graph.has_node(name_1) or not full_graph.has_node(name_2):
                continue
            if full_graph.has_edge(name_1, name_2):
                if full_graph[name_1][name_2]['metric'] > metric:
                    nx.set_edge_attributes(full_graph,{(name_1,name_2):metric}, 'metric')
            else:
                full_graph.add_edge(name_1, name_2, metric=metric)  




def _convert_to_dict(molecules: Union[List[str], np.ndarray, pd.core.series.Series, Dict[str,str]],
                    labels: Union[List[str], np.ndarray, pd.core.series.Series, Dict[str,str]] = None,
                    priority: Union[List[str], np.ndarray, pd.core.series.Series, Dict[str,str]] = None,
                    ) -> Tuple[Dict[str,str], Union[None, Dict[str,str]], Union[None, Dict[str,str]]]:
    '''
    For simplicity, we process all input data as dicts internally. Do not allow user to mix
    input types between dicts and arrays/series.
    '''

    # ensure that there is no dict/list mixup in the inputs. Could handle interally,
    # but probably better to force the user to avoid spurious mismatch errors.
    if type(molecules) == dict:
        if labels is not None and type(labels) != dict:
            raise ValueError('When providing molecules as dictionary, also provide labels as dictionary.')
        if priority is not None and type(priority) != dict:
            raise ValueError('When providing molecules as dictionary, also provide priorities as dictionary.')

        return molecules, labels, priority
    
    if type(molecules) == pd.core.series.Series:
        if labels is not None and type(labels) !=  pd.core.series.Series:
            raise ValueError('When providing molecules as series, also provide labels as series.')
        if priority is not None and type(priority) !=  pd.core.series.Series:
            raise ValueError('When providing molecules as series, also provide priorities as series.')

        labels = labels.to_dict() if labels is not None else None
        priority = priority.to_dict() if priority is not None else None
        return molecules.to_dict(), labels, priority

    elif type(molecules) in [np.ndarray, pd.core.series.Series, list]:
        if labels is not None and type(labels) not in [np.ndarray,  list]:
            raise ValueError('When providing molecules as array, also provide labels as array.')
        if priority is not None and  type(priority) not in [np.ndarray, list]:
            raise ValueError('When providing molecules as array, also provide priorities as array.')

        mol_dict = {}
        lab_dict = {} if labels is not None else None
        pri_dict = {} if priority is not None else None
        for i in range(len(molecules)):
            mol_dict[f'seq_{i}'] = molecules[i]
            if labels is not None:
                lab_dict[f'seq_{i}'] = labels[i]
            if priority is not None:
                pri_dict[f'seq_{i}'] = priority[i]

        return mol_dict, lab_dict, pri_dict


def train_test_validation_split(molecules: Union[List[str], np.ndarray, Dict[str,str]], 
                     labels: Union[List[str], np.ndarray, Dict[str,str]] = None,
                     priority: Union[List[str], np.ndarray, Dict[str,str]] = None,
                     test_size: float = 0.15, 
                     valid_size: float = 0,
                     threshold: float = 0.3,
                     initialization_mode: str = 'slow-nn',
                     no_moving: bool = False,
                     remove_same: bool = False,
                     save_checkpoint_path: str = None,
                     triangular: bool = False,
                     edge_file: str = None,
                     metric_column: str = None,
                     verbose: bool = False
                     ) -> List[Iterable]:

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ModuleNotFoundError:
        raise ImportError("This function requires RDKit to be installed.")

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

    original_type = type(molecules)
    molecules, labels, priority = _convert_to_dict(molecules, labels, priority)

    # make the graph
    full_graph, part_graph, labels = load_entities(molecules, labels, priority)
    for l in labels:
        """ Find the expected number of entities labelled l in any partition """
        labels[l]['lim'] = labels[l]['num']//partitions

    threshold = 1- threshold
    # add the edges
    compute_fingerprint_tanimoto_distances(full_graph, molecules, threshold)
    print("Full graph nr. of edges:", full_graph.number_of_edges())


    # run graph-part
    config = {
        "threshold": threshold,
        "partitions": partitions,
        "out_file": "graphpart_python", # not used anyway
        "priority_name": "priority" if priority is not None else None,
        "labels_name": "label" if labels is not None else None,
        "initialization_mode": initialization_mode,
        "no_moving": no_moving,
        "remove_same": remove_same,
        "test_ratio": test_size,
        "val_ratio": valid_size,
        "save_checkpoint_path": save_checkpoint_path,
        "triangular": triangular,
        "edge_file": edge_file,
        "metric_column": metric_column,
        "allow_moving": not no_moving, # silly conversions because in the CLI we want to have those default-false.
        "removal_type": not remove_same,
    }

    partition_assignment_df = partition_and_remove(full_graph, part_graph, labels, json_dict={}, threshold=threshold, config=config, verbose=verbose)

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


def stratified_k_fold(molecules: Union[List[str], np.ndarray, Dict[str,str]], 
                     labels: Union[List[str], np.ndarray, Dict[str,str]] = None,
                     priority: Union[List[str], np.ndarray, Dict[str,str]] = None,
                     partitions: int = 5,
                     threshold: float = 0.3,
                     initialization_mode: str = 'slow-nn',
                     no_moving: bool = False,
                     remove_same: bool = False,
                     save_checkpoint_path: str = None,
                     triangular: bool = False,
                     edge_file: str = None,
                     metric_column: str = None,
                     verbose: bool = False
                     ) -> List[Iterable]:

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ModuleNotFoundError:
        raise ImportError("This function requires RDKit to be installed.")


    original_type = type(molecules)
    molecules, labels, priority = _convert_to_dict(molecules, labels, priority)

    # make the graph
    full_graph, part_graph, labels = load_entities(molecules, labels, priority)
    for l in labels:
        """ Find the expected number of entities labelled l in any partition """
        labels[l]['lim'] = labels[l]['num']//partitions

    threshold = 1- threshold
    # add the edges
    compute_fingerprint_tanimoto_distances(full_graph, molecules, threshold)
    print("Full graph nr. of edges:", full_graph.number_of_edges())


    # run graph-part
    config = {
        "threshold": threshold,
        "partitions": partitions,
        "out_file": "graphpart_python", # not used anyway
        "priority_name": "priority" if priority is not None else None,
        "labels_name": "label" if labels is not None else None,
        "initialization_mode": initialization_mode,
        "no_moving": no_moving,
        "remove_same": remove_same,
        "test_ratio": 0,
        "val_ratio": 0,
        "save_checkpoint_path": save_checkpoint_path,
        "triangular": triangular,
        "edge_file": edge_file,
        "metric_column": metric_column,
        "allow_moving": not no_moving, # silly conversions because in the CLI we want to have those default-false.
        "removal_type": not remove_same,
    }

    partition_assignment_df = partition_and_remove(full_graph, part_graph, labels, json_dict={}, threshold=threshold, config=config, verbose=verbose)

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