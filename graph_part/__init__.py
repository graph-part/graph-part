__version__ ="1.0"

import time
from .graph_part import main

def run_graph_part():
    s = time.perf_counter()
    main()
    elapsed = time.perf_counter() - s
    print(f"{__file__} executed in {elapsed:0.2f} seconds.")