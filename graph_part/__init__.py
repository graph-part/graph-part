__version__ ="1.0"

from .cli import main
from .api import train_test_validation_split, stratified_k_fold

def run_graph_part():
    main()