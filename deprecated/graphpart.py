# -*- coding: utf-8 -*-

import sys
from datetime import datetime
from os import listdir, path
import pandas as pd 
import numpy as np
import networkx as nx
from typing import Dict, List, Tuple, Iterator
from os.path import join as path_join
from os.path import isfile, isdir
import time
from collections import Counter
from itertools import product

"""
This program partitions an entity set according to a single pairwise distance metric
and some desired threshold the partitions should fullfill.
The program seeks maximum distance between partitions. If the metric and desired 
threshold imply minimization you should designate the ('-tf', '--transformation')
parameter and specify either the one-minus or inverse transformation. The program 
converts the designated threshold the same way as the metric.
If for example a pair of entities have no distance between them (are identical)
when the metric is 1.0, and the desired threshold is 0.3, you would specify the threshold
-th 0.3 or --threshold 0.3, and then specify the desired transformation 
('-tf', '--transformation') such as -tf one-minus.
The entity/meta file can contain a priority designator and a label designator 
For example:
>entity_identifier1|exp=1|class=label1
>entity_identifier2|exp=0|class=label2
the priority name and label name should be designated using the appropriate parameters.
In this case you would specify the priority name -pn experimental and the label name 
-ln class
The following parameters are required:
    ('-mf', '--meta-file'):           Path to file with entity identifiers and 
                                      supplementary information such as labelling.
                                      Currently the interleaved fasta file format, 
                                      with | or : header separators are supported.
                                      The - header separator is untested. 
    ('-ef', '--edge-file'):           Path to a comma separated file containing 
                                      pairwise metrics, the first two columns should 
                                      contain entity identifiers specified in the 
                                      --meta-file.
    ('-th', '--threshold'):           The desired threshold, should be within the
                                      bounds defined by the metric
    ('-of', '--out-file'):            The path you want to write the partitioning to
The following parameters are technically optional, but most likely desired:
    ('-mc', '--metric-column'):       The 1-indexed or one-based indexing number,
                                      from left-to-right, specifying which column
                                      in --edge-file contains the desired metric.
                                      Left unspecified this is assumed to be 3.
    ('-tf', '--transformation'):      If some transformation is required or desired
                                      you can specify one of the following:
                                      one-minus, inverse, square or log.
    ('-pa', '--partitions'):          The number of partitions to create.
                                      By default this is set to 5.
    ('-pn', '--priority-name'):       The name of the priority in the meta file
    ('-ln', '--labels-name'):         The name of the label in the meta file
The following parameters are optional, the non default values are experimental:
    ('-im', '--initialization-mode'): 0 if you want to use the slow restricted nearest
                                      neighbour linkage, 1 if you want to use the
                                      fast restricted nearest neighbour linkage or 
                                      2 if you do not want any initialization.
                                      The default is 0.
    ('-am', '--allow-moving'):        1 allows the removing procedure to relocate
                                      entities if it finds more within threshold 
                                      neighbours in another partition. 0 dissallows 
                                      this behavior. 1 is the default.
    ('-rt', '--removal-type'):        1 if you want to remove based on the number
                                      of within threshold neighbours in other
                                      partitions. 0 if you also want to remove
                                      based on the within threshold interconnectivity
                                      between within threshold neighbours in other
                                      partitions and within threshold neighbours 
                                      in the same partition. The default is 1.
Example usage:
python graph_part.py 
    -mf ../../data/sequences/raw_data.fasta 
    -ef ../../data/edgefiles/ggs_default_gpia2_edges 
    -mc 3 
    -th 0.31 
    -pa 5 
    -tf one-minus 
    -pn experimental
    -ln positive 
    -of gpia_gps_part_identity2.csv
Author(s): Magnús Halldór Gíslason.
           Developed in close collaboration with 
           José Juan Almagro Armenteros and Henrik Nielsen.
For bug reporting contact mhgislason@gmail.com
(CC BY-NC 4.0) 
You are free to:
Share — copy and redistribute the material in any medium or format
Adapt — remix, transform, and build upon the material
Under the following terms:
Attribution — You must give appropriate credit, provide a link to the license, and indicate if changes were made. 
You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.
NonCommercial — You may not use the material for commercial purposes.
No additional restrictions — You may not apply legal terms or technological measures that legally restrict others 
from doing anything the license permits.
see:
https://creativecommons.org/licenses/by-nc/4.0/
https://creativecommons.org/licenses/by-nc/4.0/legalcode
"""

## Metric transformations should be defined here
TRANSFORMATIONS = {
    'one-minus': lambda x: 1-x, 
    'inverse': lambda x: 1/x if x > 0 else float('Inf'), 
    'square': lambda x: x**2,
    'log': lambda x: np.log(x),
    'none': lambda x: x,
    'None': lambda x: x,
    None: lambda x: x
}

## Initialization modes should be defined here, the selected mode is submitted as an integer value
INITIALIZATION_MODES = {
    0: 'slow-nn', # Full restricted nearest neighbour linkage
    1: 'fast-nn', # Partial restricted nearest neighbour linkage
    2: 'simple' # No initialization, all samples are regarded as a single 
}



def process_csv(line: str) -> Tuple[str, dict]:
    """ NOT IMPLEMENTED """
    yield None, None

def process_fasta(line: str, priority_name:str, labels_name:str, labels: dict) -> Tuple[str, dict]:
    """ Processes a fasta header lines or fasta meta data lines, if you will.
        Only supports interleaved fasta with > initialized header lines.
        Separate metadata with pipes | or colons :, whitespace between separators are ignored. 
        Some fasta files use the dash - as a separator. The dash is also used as an isoform indicator
        in accension numbers the dash separator. The option space dash space or [ - ] is implemented, but untested. """
    spl = line.strip().split('|')
    if ' - ' in spl[0]:
       spl = line.strip().split(' - ')
    if ':' in spl[0]:
        spl = line.strip().split(':')

    AC = spl[0].strip()[1:]
    priority = False
    label = '0'
    for s in spl[1:]:
        if '=' in s:
            param_spl = s.split('=')
            if param_spl[0] == priority_name:
                try:
                    priority = int(param_spl[1])==1
                except ValueError or TypeError:
                    raise TypeError("The input interpreted as priority designation did not conform as expected. Value interpeted: %r, Line: %r" % (param_spl[1], line))
            elif param_spl[0] == labels_name:
                label = str(param_spl[1].strip())

    if label not in labels:
        labels[label] = {'val':len(labels), 'num':0}
    labels[label]['num'] += 1

    label = labels[label]['val']

    node_data = {
        'priority': priority,
        'label-val': label
    }
    return AC, node_data


def load_entities(entity_fp: str, priority_name: str, labels_name: str):
    part_graph = nx.Graph()
    full_graph = nx.Graph()

    labels = {}
    with open(entity_fp) as inf:
        processing_as = None
        for line in inf:
            if '>' in line and processing_as != 'csv':
                AC, node_data = process_fasta(line, priority_name, labels_name, labels)
                processing_as = 'fasta'
            elif processing_as == 'fasta':
                continue
            else:
                AC, node_data = process_csv(line)
                processing_as = 'csv'

            full_graph.add_node(AC)
            nx.set_node_attributes(full_graph, {AC:node_data})
            
            part_graph.add_node(AC)

    return full_graph, part_graph, labels


def load_edges(edge_fp: str, 
               full_graph: nx.classes.graph.Graph, 
               part_graph: nx.classes.graph.Graph,
               tranformation: str,
               threshold: float,
               metric_column: int):
    with open(edge_fp) as inf:
        #pdb.set_trace()
        for line_nr, line in enumerate(inf):
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
             



def partition_assignment(cluster_vector, label_vector, n_partitions, n_class):
    ''' Function to separate proteins into N partitions with balanced classes 
        Courtesy of José Juan Almagro Armenteros '''
    
    # Unique cluster number
    u_cluster = np.unique(cluster_vector)
    
    # Initialize matrices
    loc_number = np.ones((n_partitions,n_class))
    cl_number = np.zeros(cluster_vector.shape[0])
    
    for i in u_cluster:
        # Extract the labels for the proteins in that cluster
        positions = np.where(cluster_vector == i)
        cl_labels = label_vector[positions]
        
        # Count number of each class
        u, count = np.unique(cl_labels, return_counts=True)
        
        u = u.astype(np.int32)
        temp_loc_number = np.copy(loc_number)
        temp_loc_number[:,u] += count
        loc_per = loc_number/temp_loc_number
        best_group = np.argmin(np.sum(loc_per,axis=1))
        loc_number[best_group,u] += count
        
        # Store the selected partition
        cl_number[positions] = best_group
    
    return cl_number
        

def partition_data(full_graph: nx.classes.graph.Graph, 
                   part_graph: nx.classes.graph.Graph,
                   labels: dict,
                   threshold: float,
                   nr_of_parts: int,
                   mode: int):
    part_size = full_graph.number_of_nodes()//nr_of_parts

    label_limits = np.array([x[1]['lim'] for x in sorted(labels.items(), key=lambda x:x[1]['val'] )])
    print(part_size, label_limits)
    
    count = 0
    ## AC, cluster, label
    acs = []
    clusters = []
    labels = []
    print("Initialization mode", INITIALIZATION_MODES[mode])

    ## Initialize the initialization
    for ind, AC in enumerate(part_graph.nodes()):
        label_counts = np.zeros(len(label_limits), dtype=int)
        label_counts[full_graph.nodes[AC]['label-val']] = 1
        cluster_nr = ind
        
        if mode == 2:
            d = full_graph.nodes[AC]
            acs.append(AC)
            clusters.append(cluster_nr)
            labels.append(d['label-val'])
        nx.set_node_attributes(part_graph, {
            AC:{
                'cluster': cluster_nr,
                'C-size': 1,
                'label-counts': label_counts
            }
        })

    ## Restricted closest neighbour linkage
    if mode < 2:
        ## Linking entities, if restrictions allow
        for qry, lib, data in sorted(full_graph.edges(data=True), key=lambda x: x[2]['metric']):
            if data['metric'] > threshold:
                ## No need to continue if threshold reached.
                break
            
            if part_graph.has_edge(qry, lib):
                ## Update edge if it exists
                if part_graph[qry][lib]['metric'] > data['metric']:
                    nx.set_edge_attributes(part_graph,{(qry,lib):data['metric']}, 'metric')
                continue
            
            if part_graph.nodes[qry]['cluster'] == part_graph.nodes[lib]['cluster']:
               continue

            ## RESTRICTIONS!
            if part_graph.nodes[qry]['C-size'] >= part_size:
                continue
            if part_graph.nodes[lib]['C-size'] >= part_size:
                continue
            if (part_graph.nodes[qry]['label-counts']+part_graph.nodes[lib]['label-counts'] >= label_limits).any():
                continue
            
            ## Add edge and update mini-cluster
            attr = part_graph.nodes[qry]
            if attr['cluster'] != part_graph.nodes[lib]['cluster'] or INITIALIZATION_MODES[mode] == 'fast-nn':
                attr['C-size'] += part_graph.nodes[lib]['C-size']
                attr['label-counts'] += part_graph.nodes[lib]['label-counts']
            nx.set_node_attributes(part_graph, {qry:attr})
            part_graph.add_edge(qry, lib, metric=data['metric'])

            for descendant in nx.descendants(part_graph, qry):
                nx.set_node_attributes(part_graph, {descendant:attr})

            count += 1
            if count % 10000 == 0:
                print("edges:", part_graph.number_of_edges()) 
                print (attr, data)
        #with open('graph_part_miniclusters.txt','w+') as outf:
        for cc_nr, cc in enumerate(nx.connected_components(part_graph)):
            for n in cc:
                d = full_graph.nodes[n]
                acs.append(n)
                clusters.append(cc_nr)
                labels.append(d['label-val'])
                    #outf.write("%s,%d\n" % (n,cc_nr))

    acs = np.array(acs)
    clusters = np.array(clusters)
    labels = np.array(labels)
    print(len(np.unique(labels)))
    partitioning = partition_assignment(clusters, labels, nr_of_parts, len(np.unique(labels)))
    for ind, p in enumerate(partitioning):
        attr = part_graph.nodes[acs[ind]]
        attr['cluster'] = p
        nx.set_node_attributes(part_graph,{acs[ind]:attr})


def remover(full_graph: nx.classes.graph.Graph, 
            part_graph: nx.classes.graph.Graph, 
            threshold:float, 
            move_to_most_neighbourly:bool = True, 
            ignore_priority:bool = True,
            simplistic_removal:bool = True):
    
    print("Min-threshold", "\t", "#Entities", "\t", "#Edges", "\t", "Connectivity", "\t", "#Problematics", "\t", "#Relocated", "\t", "#To-be-removed")
    removing_round = 0
    while True:
        between_connectivity = {}
        min_oc_wth= 1
        number_moved = 0
        for n,d in full_graph.nodes(data=True):
            neighbours = nx.neighbors(full_graph,n)
            neighbour_clusters = Counter((part_graph.nodes[nb]['cluster'] for nb in nx.neighbors(full_graph,n) if full_graph[n][nb]['metric'] < threshold))
            cluster = part_graph.nodes[n]['cluster']
            nb_sc_wth = []
            nb_oc_wth = []
            
            ## FIRST MOVE NODE TO CLUSTER WITH MOST NEIGHBOURS
            if move_to_most_neighbourly:
                if len(neighbour_clusters) > 0:
                    most_neighbourly_cluster = max(neighbour_clusters.items(), key=lambda x:x[1])[0]
                    if most_neighbourly_cluster != cluster:
                        part_graph.nodes[n]['cluster'] = most_neighbourly_cluster
                        number_moved += 1
            
            if ignore_priority and full_graph.nodes[n]['priority']:
                between_connectivity[n] = 0
                continue

            for neighbour in neighbours:
                nb_cluster = part_graph.nodes[neighbour]['cluster']
                if nb_cluster == cluster and full_graph[n][neighbour]['metric'] < threshold and not simplistic_removal:
                    ## The more complex removal criterion. This was worse for GPI anchors. Haven't checked for Subnuclear
                    nb_sc_wth.append(full_graph[n][neighbour]['metric'])
                elif nb_cluster != cluster and full_graph[n][neighbour]['metric'] < threshold:
                    min_oc_wth = min(min_oc_wth, full_graph[n][neighbour]['metric'])
                    nb_oc_wth.append(full_graph[n][neighbour]['metric'])

            ## The more complex additional connectivity criterion, will be empty with simple_removal == True        
            distances = [x for x in (sum(x_) for x_ in product(nb_sc_wth, nb_oc_wth)) if x >= threshold] 
            
            distances += nb_oc_wth

            between_connectivity[n] = len(distances)
            
        nx.set_node_attributes(full_graph, between_connectivity, 'between_connectivity')
        bc_sum = np.sum(np.fromiter((d['between_connectivity'] for n,d in full_graph.nodes(data=True)),int))
        bc_count = np.sum(np.fromiter((1 for n,d in full_graph.nodes(data=True) if d['between_connectivity'] > 0),int))

        removing_round += 1
        number_to_remove = int(bc_count*np.log10(removing_round)/100)+1 # int(bc_count*0.01)+1
        ## Remove 1% + 1 of the most problematic entities
        remove_these = [x[0] for x in sorted(((n,d['between_connectivity']) for n,d in full_graph.nodes(data=True) if d['between_connectivity'] > 0), key=lambda x:x[1], reverse=True)[:number_to_remove]]
        
        print(round(min_oc_wth,7), "\t\t", full_graph.number_of_nodes(), "\t\t", full_graph.number_of_edges(), "\t\t", bc_sum, "\t\t", bc_count, "\t\t", number_moved, "\t\t", len(remove_these))

        full_graph.remove_nodes_from(remove_these)
        # If we've removed the last problematic entities, we stop
        if full_graph.number_of_nodes()==0 or bc_sum==0 or len(remove_these) == bc_count:
            break

def score_partitioning(df:pd.core.frame.DataFrame) -> float:
    s0 = df.shape[0]
    s1 = df.shape[1]
    return float((df.product(axis=1)**(1/s1)).product()**(1/s0))

def display_results(
    part_graph: nx.classes.graph.Graph, 
    full_graph: nx.classes.graph.Graph,
    labels: dict,
    nr_of_parts: int) -> Tuple[pd.core.frame.DataFrame, pd.core.frame.DataFrame]:
    """ """
    df = pd.DataFrame(((d) for n,d in full_graph.nodes(data=True)))
    df['cluster'] = [part_graph.nodes[n]['cluster'] for n in full_graph.nodes()]

    # It can happen that removal completely removed one partition.
    # In this case, we need to report back an error
    if len(df['cluster'].unique()) < nr_of_parts:
        error_string = f'''
        Impossible to generate the desired {nr_of_parts} partitions at the current partitioning threshold.
        Removal of sequences to achieve separation results in loss of {nr_of_parts-len(df['cluster'].unique())} complete partitions.
        '''
        raise RuntimeError(error_string)

    df['AC'] = [n for n in full_graph.nodes()]
    result = df.groupby(['cluster','label-val'])['AC'].count().reset_index().pivot_table(values='AC',columns=['label-val'],index=['cluster']).T
    df.set_index('AC', inplace=True)
    result['label'] = ''
    for l in labels:
        result.loc[labels[l]['val'], 'label'] = l
        result['mean'] = result[list(range(nr_of_parts))].mean(axis=1)
        result['count'] = result[list(range(nr_of_parts))].sum(axis=1)
    print(result)
    print()
    print("Partitioning score:", score_partitioning(result[range(nr_of_parts)]))
    print()
    return df, result
    

def removal_needed(
    part_graph: nx.classes.graph.Graph, 
    full_graph: nx.classes.graph.Graph,
    threshold: float) -> bool:
    """ """
    min_between = float('Inf')
    for qry, lib, fed in full_graph.edges(data=True):
        if part_graph.nodes[qry]['cluster'] != part_graph.nodes[lib]['cluster']:
            if fed['metric'] < min_between:
                min_between = fed['metric']
            if min_between < threshold:
                print ("! ", qry, lib, fed, " !")
                return True
    return False

def main():
    
    entity_fp = None
    edge_fp = None
    threshold = None
    nr_of_parts = 5
    transformation = None
    priority_name = None
    labels_name = None
    outfile = None
    allow_moving = True
    removal_type = True
    mode = 0
    metric_column = 3
    if len(sys.argv) > 1:
        args = (x for x in sys.argv[1:])
        for arg in args:
            parts = arg
            if ':' in arg:
                parts = arg.split(':')
            elif '=' in arg:
                parts = arg.split('=')
            elif arg in ('-mf', '--meta-file'):
                entity_fp = next(args, None)
            elif arg in ('-ef', '--edge-file'):
                edge_fp = next(args, None)
            elif arg in ('-mc', '--metric-column'):
                metric_column = next(args, None)
            elif arg in ('-th', '--threshold'):
                threshold = next(args, None)
            elif arg in ('-pa', '--partitions'):
                nr_of_parts = next(args, None)
            elif arg in ('-tf', '--transformation'):
                transformation = next(args, None)
            elif arg in ('-pn', '--priority-name'):
                priority_name = next(args, None)
            elif arg in ('-ln', '--labels-name'):
                labels_name = next(args, None)
            elif arg in ('-of', '--out-file'):
                outfile = next(args, None)
            elif arg in ('-im', '--initialization-mode'):
                mode = next(args, None)
            elif arg in ('-am', '--allow-moving'):
                allow_moving = next(args, None)
            elif arg in ('-rt', '--removal-type'):
                removal_type = next(args, None)
            
            if len(parts) == 2:
                if parts[0] in ('-mf', '--meta-file'):
                    entity_fp = parts[1]
                elif parts[0] in ('-ef', '--edge-file'):
                    edge_fp = parts[1]
                elif parts[0] in ('-mc', '--metric-column'):
                    metric_column = parts[1]
                elif parts[0] in ('-th', '--threshold'):
                    threshold = parts[1]
                elif parts[0] in ('-pa', '--partitions'):
                    nr_of_parts = parts[1]
                elif parts[0] in ('-tf', '--transformation'):
                    transformation = parts[1]
                elif parts[0] in ('-pn', '--priority-name'):
                    priority_name = parts[1]
                elif parts[0] in ('-ln', '--labels-name'):
                    labels_name = parts[1]
                elif parts[0] in ('-of', '--out-file'):
                    outfile = parts[1]
                elif parts[0] in ('-im', '--initialization-mode'):
                    mode = parts[1]
                elif parts[0] in ('-am', '--allow-moving'):
                    allow_moving = parts[1]
                elif parts[0] in ('-rt', '--removal-type'):
                    removal_type = parts[1]
    
    ## Validate Meta file path
    if type(entity_fp) != str:
        raise TypeError("Meta file path (-mf/--meta-file) not provided, improper or nonexistent.") 
    elif not isfile(entity_fp):
        raise ValueError("Meta file path (-mf/--meta-file) improper or nonexistent.") 
    
    ## Validate Edge file path
    if type(edge_fp) != str:
        raise TypeError("Edge file path (-ef/--edge-file) not provided, improper or nonexistent.") 
    elif not isfile(edge_fp):
        raise ValueError("Edge file path (-ef/--edge-file) improper or nonexistent.") 
    
    ## Validate out file path
    if type(outfile) != str:
        raise TypeError("Out file path (-of/--out-file) not provided or improper.") 
    else:
        try:
            with open(outfile, 'w+') as outf:
                pass
        except:
            raise ValueError("Entity file path (-of/--out-file) improper or nonexistent.") 
    
    ## Validate transformation
    if transformation is not None and transformation not in TRANSFORMATIONS:
        raise ValueError("The transformation entered (-tf/--transformation) %s, is not recognized. Please select one of the following: %s." % (transformation, ', '.join(list(TRANSFORMATIONS.keys())[:-2])))
    
    ## Validate mode
    try:
        mode = int(mode)
    except ValueError or TypeError:    
        raise ValueError("The initialization mode entered (-im/--initialziation-mode) %s, is not recognized. Please select one of the following: %s." % (mode, ', '.join([str(x) for x in INITIALIZATION_MODES.keys()])))
    if mode not in INITIALIZATION_MODES:
        raise ValueError("The initialization mode entered (-im/--initialziation-mode) %s, is not recognized. Please select one of the following: %s." % (mode, ', '.join([str(x) for x in INITIALIZATION_MODES.keys()])))
    
    ## Validate threshold
    try:
        threshold = TRANSFORMATIONS[transformation](float(threshold))
    except ValueError or TypeError:
        raise TypeError("Failed to interpret the entered threshold (-th/--threshold)%r. Please enter a value interpretable as a real value." % (threshold))
        return
    
    ## Validate number of partitions
    try:
        nr_of_parts = int(nr_of_parts)
    except ValueError or TypeError:
        raise TypeError("Failed to interpret the entered number of partitions  (-pa/--partitions) %r. Please enter a (reasonable) integer value." % (nr_of_parts))
        
    ## Validate allow moving 
    try:
        allow_moving = int(allow_moving) == 1
    except ValueError or TypeError:
        raise TypeError("Failed to interpret allow moving  (-am/--allow-moving) %r. Please enter 0 to dissallow moving otherwise don't specify or enter 1" % (allow_moving))
        
    ## Validate removal type
    try:
        removal_type = int(removal_type) == 1
    except ValueError or TypeError:
        raise TypeError("Failed to interpret removal type  (-rt/--removal-type) %r. Please enter 0 for alternative removal otherwise don't specify or enter 1" % (removal_type))
    
    ## Validate Edge file metric column
    try:
        metric_column = int(metric_column)-1
        print('metric column', metric_column)
    except ValueError or TypeError:
        raise TypeError("Failed to interpret the entered edge file metric column (-mc/--metric-column) %r. Please enter an integer value representing the 1-based, left-to-right, column index for the desired metric in the supplied edge list file." % (nr_of_parts))

    ## End of parameter validation

    ## Processing starts here:

    ## Load entities/samples as networkx graphs. labels contains label metadata.
    full_graph, part_graph, labels = load_entities(entity_fp, priority_name, labels_name)

    for l in labels:
        """ Find the expected number of entities labelled l in any partition """
        labels[l]['lim'] = labels[l]['num']//nr_of_parts

    ## Let's see the initial label distribution
    print(pd.DataFrame(labels).T)

    ## Load the edges
    load_edges(edge_fp, full_graph, part_graph, transformation, threshold, metric_column)

    ## Let's look at the number of edges
    print("Full graph nr. of edges:", full_graph.number_of_edges())
    
    ## Finally, let's partition this
    partition_data(full_graph, part_graph, labels, threshold, nr_of_parts, mode)

    df, _ = display_results(part_graph, full_graph, labels, nr_of_parts)
    df.to_csv(outfile + "pre-removal")
    print('Currently have this many samples:', full_graph.number_of_nodes())
    
    ## Check if we need to remove any
    if removal_needed(part_graph, full_graph, threshold):     
        print('Need to remove! Currently have this many samples:', full_graph.number_of_nodes())

        remover(full_graph, part_graph, threshold, allow_moving, True, removal_type)    

    if removal_needed(part_graph, full_graph, threshold):   
        print('Need to remove priority! Currently have this many samples:', full_graph.number_of_nodes())
        remover(full_graph, part_graph, threshold, allow_moving, False, removal_type)    

    print('After removal we have this many samples:', full_graph.number_of_nodes())


    df, _ = display_results(part_graph, full_graph, labels, nr_of_parts)

    if removal_needed(part_graph, full_graph, threshold):
        print ("Something is wrong! Removal still needed!")

    ## clustering to outfile. This will probably change...
    df.to_csv(outfile)

if __name__ == "__main__":
    s = time.perf_counter()
    main()
    elapsed = time.perf_counter() - s
    print(f"{__file__} executed in {elapsed:0.2f} seconds.")