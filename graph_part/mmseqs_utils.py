import subprocess
import os
import shutil
from .transformations import TRANSFORMATIONS
import networkx as nx
from tqdm.auto import tqdm

# mmseqs createdb data/netgpi_dataset.fasta temp/netgpi_db
# mmseqs prefilter -s 7.5 temp/netgpi_db temp/netgpi_db temp/netgpi_pref
# mmseqs align temp/netgpi_db temp/netgpi_db temp/netgpi_pref temp/netgpi_align_db
# mmseqs convertalis temp/netgpi_db temp/netgpi_db temp/netgpi_align_db temp/alignments.tab # column 3 contains sequence identity.
# paste -d "," <(cut -f1 -d "|" temp/alignments.tab) <(cut -f2 temp/alignments.tab | cut -f1 -d "|") <(cut -f3 temp/alignments.tab) >netgpi_mmseqs_edgelist.csv



def generate_edges_mmseqs(entity_fp: str, 
                  full_graph: nx.classes.graph.Graph, 
                  tranformation: str,
                  threshold_transformed: float,
                  threshold_original: float = None, # use original threshold (before one-minus) to pass to mmseqs
                  denominator: str = 'longest',
                  delimiter: str = '|',
                  is_nucleotide: bool = False,
                  use_prefilter: bool = False,
                  ) -> None:


    if shutil.which('mmseqs') is None:
        print('MMseqs2 was not found. Please run `conda install -c conda-forge -c bioconda mmseqs2`')
        exit()

    os.makedirs('temp', exist_ok=True)

    # Run all mmseqs ops to get a tab file that contains the alignments.
    typ = '2' if is_nucleotide else '1'
    subprocess.run(['mmseqs', 'createdb', '--dbtype', typ, entity_fp, 'temp/seq_db'])

    # However, this function will not work with nucleotidenucleotide searches, 
    # since we need to have a valid diagonal for the banded alignment.
    if is_nucleotide or use_prefilter:
        subprocess.run(['mmseqs', 'prefilter', '-s', '7.5', 'temp/seq_db', 'temp/seq_db', 'temp/pref'])
    else:
        subprocess.run(['mmseqs_fake_prefilter.sh', 'temp/seq_db', 'temp/seq_db', 'temp/pref', 'seq_db'])

    # 0: alignment length 1: shorter, 2: longer sequence
    id_mode = {'n_aligned':'0', 'shortest':'1', 'longest':'2'}[denominator]
    
    command = ['mmseqs', 'align',  'temp/seq_db', 'temp/seq_db', 'temp/pref', 'temp/align_db', '--alignment-mode', '3', '-e', 'inf', '--seq-id-mode', id_mode]
    if threshold_original is not None:
        command = command + ['--min-seq-id', str(threshold_original)]
    subprocess.run(command)
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
            if metric > threshold_transformed:
                continue
            if not full_graph.has_node(this_qry) or not full_graph.has_node(this_lib):
                continue
            if full_graph.has_edge(this_qry, this_lib):
                if full_graph[this_qry][this_lib]['metric'] > metric:
                    #nx.set_edge_attributes(full_graph,{(this_qry,this_lib):metric}, 'metric')
                    full_graph.add_edge(this_qry, this_lib, metric=metric) #Notes: Adding an edge that already exists updates the edge data. 
            else:
                full_graph.add_edge(this_qry, this_lib, metric=metric)  

    shutil.rmtree('temp')
