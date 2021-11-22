'''
Parsing functions for precomputed similarities.
'''
import networkx as nx
from .transformations import TRANSFORMATIONS
from tqdm import tqdm

def load_edge_list(edge_fp: str, 
               full_graph: nx.classes.graph.Graph, 
               tranformation: str,
               threshold: float,
               metric_column: int):
    '''
    Load edges form a precomputed edge list saved as .csv
    Expects the names of the nodes in columns 0 and 1, the
    metric in metric_column.
    '''
    with open(edge_fp) as inf:
        #pdb.set_trace()
        for line_nr, line in tqdm(enumerate(inf)):
            spl = line.strip().split(',')
            if len(spl) < 3:
                raise ValueError("""
                Edge list file does not contain at least three comma 
                separated columns. The first two columns should contain
                entity identifiers and the third should contain the
                metric to partition by.
                """)
            
            this_qry = spl[0]
            this_lib = spl[1]
            try:
                metric = TRANSFORMATIONS[tranformation](float(spl[metric_column]))
            except ValueError or TypeError:
                raise TypeError("Failed to interpret the metric column value %r. Please ensure that the edge list file is correctly formatted and that the correct column is specified." % (spl[1]))
            
            if this_qry == this_lib:
                continue
            if metric > threshold:
                continue
            if not full_graph.has_node(this_qry) or not full_graph.has_node(this_lib):
                continue
            if full_graph.has_edge(this_qry, this_lib):
                if full_graph[this_qry][this_lib]['metric'] > metric:
                    nx.set_edge_attributes(full_graph,{(this_qry,this_lib):metric}, 'metric')
            else:
                full_graph.add_edge(this_qry, this_lib, metric=metric)  