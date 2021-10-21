import subprocess
import os
import shutil
from .transformations import TRANSFORMATIONS
import networkx as nx
from tqdm import tqdm

# mmseqs createdb data/netgpi_dataset.fasta temp/netgpi_db
# mmseqs prefilter -s 7.5 temp/netgpi_db temp/netgpi_db temp/netgpi_pref
# mmseqs align temp/netgpi_db temp/netgpi_db temp/netgpi_pref temp/netgpi_align_db
# mmseqs convertalis temp/netgpi_db temp/netgpi_db temp/netgpi_align_db temp/alignments.tab # column 3 contains sequence identity.
# paste -d "," <(cut -f1 -d "|" temp/alignments.tab) <(cut -f2 temp/alignments.tab | cut -f1 -d "|") <(cut -f3 temp/alignments.tab) >netgpi_mmseqs_edgelist.csv



def generate_edges_mmseqs(entity_fp: str, 
                  full_graph: nx.classes.graph.Graph, 
                  tranformation: str,
                  threshold: float,
                  delimiter: str = '|',
                  is_nucleotide: bool = False,
                  ) -> None:

    os.makedirs('temp', exist_ok=True)

    # Run all mmseqs ops to get a tab file that contains the alignments.
    typ = '2' if is_nucleotide else '1'
    subprocess.run(['mmseqs', 'createdb', '--dbtype', typ, entity_fp, 'temp/seq_db'])
    subprocess.run(['mmseqs', 'prefilter', '-s', '7.5', 'temp/seq_db', 'temp/seq_db', 'temp/pref'])
    subprocess.run(['mmseqs', 'align', 'temp/seq_db', 'temp/seq_db', 'temp/pref', 'temp/align_db'])
    subprocess.run(['mmseqs', 'convertalis', 'temp/seq_db', 'temp/seq_db', 'temp/align_db', 'temp/alignments.tab'])

    # Read the result
    with open('temp/alignments.tab') as inf:
        for line_nr, line in tqdm(enumerate(inf)):
            spl = line.strip().split('\t')

            this_qry = spl[0].split(delimiter)[0]
            this_lib = spl[1].split(delimiter)[0]
            ident = float(spl[2])

            try:
                metric = TRANSFORMATIONS[tranformation](ident)
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

    shutil.rmtree('temp')
