## Metric transformations should be defined here
import numpy as np

TRANSFORMATIONS = {
    'one-minus': lambda x: 1-x, 
    'inverse': lambda x: 1/x if x > 0 else float('Inf'), 
    'square': lambda x: x**2,
    'log': lambda x: np.log(x),
    'none': lambda x: x,
    'None': lambda x: x,
    None: lambda x: x
}


INVERSE_TRANSFORMATIONS = {
    'one-minus': lambda x: 1-x, 
    'inverse': lambda x: 1/x if x > 0 else float('Inf'), 
    'square': lambda x: np.sqrt(x), #x**2,
    'log': lambda x: np.exp(x), #np.log(x),
    'none': lambda x: x,
    'None': lambda x: x,
    None: lambda x: x
}
