'''
Command line interface for Graph-Part.
'''
import argparse
import os
from .transformations import TRANSFORMATIONS
from .train_val_test_split import check_train_val_test_args
from .graph_part import run_partitioning

#TODO check all help strings and update if needed
def get_args() -> argparse.Namespace:
    '''
    Create an ArgumentParser to process all command line arguments of Graph-Part.
    Returns the parsed arguments.
    Also runs some checks to validate argument combinations.
    '''

    parser = argparse.ArgumentParser('Graph-Part')#, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # 1. Arguments that are always required.
    # Input processing and graphpart core parameters.
    core_parser = argparse.ArgumentParser(add_help=False)
    core_parser.add_argument("-ff","--fasta-file",type=str, help='''Path to file with entity identifiers and 
                                                            supplementary information such as labelling.
                                                            Currently the interleaved fasta file format, 
                                                            with | or : header separators are supported.
                                                            The - header separator is untested. ''',
                        required=True,
                        )
    core_parser.add_argument("-th","--threshold",type=float, help='''The desired threshold, should be within the
                                                              bounds defined by the metric''',
                        required=True,
                        )
    core_parser.add_argument("-pa","--partitions",type=int, help='Number of partitions to generate.', 
                        default=5,
                        )
    core_parser.add_argument("-tf","--transformation",type=str, help='Transformation to apply to the similarity/distance metric.', 
                        choices=list(TRANSFORMATIONS.keys()), 
                        default='one-minus',
                        )
    core_parser.add_argument("-of","--out-file",type=str, help='The path you want to write the partitioning to.', default='graphpart_result.csv')


    core_parser.add_argument("-pn","--priority-name",type=str, help='The name of the priority in the meta file.', 
                        default=None,
                        )
    core_parser.add_argument("-ln","--labels-name",type=str, help='The name of the label in the meta file.', 
                        default=None,
                        )
    
    core_parser.add_argument("-im","--initialization-mode",type=str, help='Use either slow or fast restricted nearest neighbor linkage or no initialization.', 
                        default='slow-nn', 
                        choices=['slow-nn', 'fast-nn', 'simple'],
                        )
    core_parser.add_argument("-nm","--no-moving",action='store_true', help='''Disallows the removing procedure to relocate
                                                                            entities if it finds more within threshold 
                                                                            neighbours in another partition.'''
                        )
    core_parser.add_argument("-rs","--remove-same",action='store_true', help='''Activate if you also want to remove
                                      based on the within threshold interconnectivity
                                      between within threshold neighbours in other
                                      partitions and within threshold neighbours 
                                      in the same partition. The default removes based on the number
                                      of within threshold neighbours in other partitions.'''
                        )
    
    # train-val-test splits.
    core_parser.add_argument("-te","--test-ratio",type=float, default=0.0, help='The fraction of the data to use for testing. Incompatible with `partitions`.')
    core_parser.add_argument("-va","--val-ratio",type=float, default=0.0,help='The fraction of the data to use for validation. Incompatible with `partitions`.')
    
    #checkpointing
    core_parser.add_argument('--save-checkpoint-path', '-sc', type=str, default=None, help='''Path to save the computed similarities above the threshold. 
                                                                                            Can be used later in the precomputed mode.'''
                        )


    # Parsers for the different run modes.
    subparsers = parser.add_subparsers(title='modes',description='Available alignment modes.', dest='alignment_mode')
    subparsers.required = True # ugly but apparently this is how you force a subparser to be specified.
    parser_precomputed =  subparsers.add_parser('precomputed', help='Use precomputed identities.', parents=[core_parser])
    parser_needle = subparsers.add_parser('needle', help='Use EMBOSS needle alignments.', parents=[core_parser])
    parser_mmseqs2 = subparsers.add_parser('mmseqs2', help='Use MMseqs2 alignments.', parents=[core_parser])

    # 2. Arguments that are only required with precomputed metrics.
    parser_precomputed.add_argument("-ef","--edge-file",type=str, help='''Path to a comma separated file containing 
                                                            pairwise metrics, the first two columns should 
                                                            contain entity identifiers specified in the 
                                                            --fasta-file.''',
                        default=None,
                        )
    parser_precomputed.add_argument("-mc","--metric-column",type=int, help='''The 0-indexed or zero-based indexing number,
                                                            from left-to-right, specifying which column
                                                            in --edge-file contains the desired metric.
                                                            Left unspecified this is assumed to be 2.''', 
                        default=2,
                        )

    # 3. Arguments that are only required with needleall.
    parser_needle.add_argument("-dn","--denominator",type=str, help='Denominator to use for sequence identity computation.', 
                        choices=['full', 'shortest', 'longest', 'mean', 'no_gaps'], 
                        default='full',
                        )
    # Flags
    parser_needle.add_argument("-nu","--nucleotide", action='store_true', help= 'Input contains nucleotide sequences (Default is proteins).')
    parser_needle.add_argument("-tr","--triangular", action='store_true', help='Only compute triangular part of full distance matrix.')

    # optimize runtime
    parser_needle.add_argument("-nt","--threads",type=int, help='Number of threads to run in parallel.', default=1)
    parser_needle.add_argument("-nc","--chunks",type=int, help='Number of chunks to split the fasta file.', default=10)
    parser_needle.add_argument("-pm", "--parallel-mode", type=str, help='Parallelization strategy to use.',
                                choices=['multithread', 'multiprocess'], 
                                default='multithread'
                                )
    
    # customize needle
    parser_needle.add_argument('--gapopen','-gapopen', type=float, default=10, help='Passed to needle. See EMBOSS documentation.')
    parser_needle.add_argument('--gapextend','-gapextend', type=float, default=0.5, help='Passed to needle. See EMBOSS documentation.')
    parser_needle.add_argument('--endweight','-endweight', action='store_true', help='Passed to needle. See EMBOSS documentation.')
    parser_needle.add_argument('--endopen','-endopen', type=float, default=10, help='Passed to needle. See EMBOSS documentation.')
    parser_needle.add_argument('--endextend','-endextend', type=float, default=0.5, help='Passed to needle. See EMBOSS documentation.')
    parser_needle.add_argument('--matrix', '--datafile','-datafile', type=str, default='EBLOSUM62', help='Passed to needle. See EMBOSS documentation.')


    # 4. Arguments that are only required with mmseqs2.
    parser_mmseqs2.add_argument("-nu","--nucleotide", action='store_true', help= 'Input contains nucleotide sequences (Default is proteins).')
    parser_mmseqs2.add_argument("-pr","--prefilter", action='store_true', help= 'Use the mmseqs2 prefiltering procedure instead of forcing all-vs-all alignments.')
    parser_mmseqs2.add_argument("-dn","--denominator",type=str, help='Denominator to use for sequence identity computation.', 
                        choices=['shortest', 'longest', 'n_aligned'], 
                        default='shortest',
                        )

    args =  parser.parse_args()


    # Perform checks
    def create_dir_or_fail(file_path: str) -> None:
        '''Create output dirs or fail if not possible. If no dir, check that current working dir is writeable or fail.'''

        out_path = os.path.dirname(file_path)
        if out_path != '': # dirname is empty when we just have a filename.
            if not os.path.exists(out_path):
                os.makedirs(out_path, exist_ok=True)
        else:
            # if we just have a filename, validate that the current working directory is writeable.
            if not os.access(os.getcwd(), os.W_OK):
                raise PermissionError(file_path)


    create_dir_or_fail(args.out_file)
    if args.save_checkpoint_path is not None:
        create_dir_or_fail(args.save_checkpoint_path)

    if args.test_ratio >0 or args.val_ratio>0:
        check_train_val_test_args(args)
        print(f'Running in train-validation-test split mode with {args.test_ratio:.0%} test and {args.val_ratio:.0%} validation.')
    
            
    return args

def main():

    
    args = get_args()
    config = vars(args)
    config['allow_moving'] = not args.no_moving
    config['removal_type'] = not args.remove_same

    run_partitioning(config, write_output_file=True, write_json_report=True)
