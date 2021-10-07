'''
Command line interface for Graph-Part.
'''
import argparse
from .transformations import TRANSFORMATIONS

#TODO check all help strings and update if needed
def get_args() -> argparse.Namespace:
    '''
    Create an ArgumentParser to process all command line arguments of Graph-Part.
    Returns the parsed arguments.
    Also runs some checks to validate argument combinations.
    '''

    parser = argparse.ArgumentParser('Graph-Part', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-ff","--fasta-file",type=str, help='''Path to file with entity identifiers and 
                                                            supplementary information such as labelling.
                                                            Currently the interleaved fasta file format, 
                                                            with | or : header separators are supported.
                                                            The - header separator is untested. ''',
                        required=True,
                        )
    parser.add_argument("-ef","--edge-file",type=str, help='''Path to a comma separated file containing 
                                                            pairwise metrics, the first two columns should 
                                                            contain entity identifiers specified in the 
                                                            --fasta-file.''',
                        default=None,
                        )
    parser.add_argument("-mc","--metric-column",type=int, help='''The 0-indexed or zero-based indexing number,
                                                            from left-to-right, specifying which column
                                                            in --edge-file contains the desired metric.
                                                            Left unspecified this is assumed to be 2.''', 
                        default=3,
                        )
    parser.add_argument("-th","--threshold",type=float, help='''The desired threshold, should be within the
                                                              bounds defined by the metric''',
                        required=True,
                        )
    parser.add_argument("-pa","--partitions",type=int, help='Number of partitions to generate.', 
                        default=5,
                        )
    parser.add_argument("-tf","--transformation",type=str, help='Transformation to apply to the similarity/distance metric.', 
                        choices=list(TRANSFORMATIONS.keys()), 
                        default=None,
                        )
    parser.add_argument("-of","--out-file",type=str, help='The path you want to write the partitioning to.', default='graphpart_result.csv')


    parser.add_argument("-pn","--priority-name",type=str, help='The name of the priority in the meta file.', 
                        default=None,
                        )
    parser.add_argument("-ln","--labels-name",type=str, help='The name of the label in the meta file.', 
                        default=None,
                        )
    
    parser.add_argument("-im","--initialization-mode",type=str, help='Use either slow or fast restricted nearest neighbor linkage or no initialization.', 
                        default='slow-nn', 
                        choices=['slow-nn', 'fast-nn', 'simple'],
                        )
    
    
    


    # Flags
    parser.add_argument("-nm","--no-moving",action='store_true', help='''Disallows the removing procedure to relocate
                                                                            entities if it finds more within threshold 
                                                                            neighbours in another partition.'''
                        )
    parser.add_argument("-rs","--remove-same",action='store_true', help='''Activate if you also want to remove
                                      based on the within threshold interconnectivity
                                      between within threshold neighbours in other
                                      partitions and within threshold neighbours 
                                      in the same partition. The default removes based on the number
                                      of within threshold neighbours in other partitions.'''
                        )

    parser.add_argument("-nt","--threads",type=int, help='Number of threads to run in parallel.', default=1)
    parser.add_argument("-nc","--chunks",type=int, help='Number of chunks to split the fasta file.', default=10)

    parser.add_argument('--load-checkpoint-path', '-lc', type=str, default=None, help='Path to save the generated graph.')
    parser.add_argument('--save-checkpoint-path', '-sc', type=str, default=None, help='Path to a previously generated graph for quickstart.')

    args =  parser.parse_args()

    ## Validate argument combinations.
    #TODO rather make those warnings
    if args.load_checkpoint_path is not None and args.threads>1:
        print('Cannot use parallel processing when starting from a precomputed graph.  --threads argument will be ignored.')

    if args.edge_file is not None and args.threads>1:
        print('Cannot use parallel processing when starting from a precomputed edge list. --threads argument will be ignored.')

    

    # prevent unintentional overwriting of previous graph. Can be very expensive to recompute for large datasets
    if args.load_checkpoint_path is not None and (args.load_checkpoint_path == args.save_checkpoint_path): 
        print(f'load-checkpoint-path and save-checkpoint-path are identical! This would overwrite the currently saved graph with a new one thresholded at {args.threshold}')
        print('Please confirm with y to continue, n to abort:')

        res = ''
        while res not in {"y", "n"}:
            res = input("(Enter y/n)").lower()
        
        if res == 'n':
            exit()
        
    
    return args


