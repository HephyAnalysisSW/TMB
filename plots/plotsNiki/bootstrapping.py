''' Implement bootstrapping (resampling) 
    https://en.wikipedia.org/wiki/Bootstrapping_(statistics)
'''

import random
import numpy as np

def boot_strapped_quantiles(  data, f, n, quantiles = [0.0445, 0.3173, 0.5, 1-0.3173, 1-0.0445]):
        ''' make n replica data sets (of the same size) by random sampling the original data with replacements.
            Evaluate a function on the replicas and compute the quantiles of the results.''' 
        replica_values = sorted(map( f,  ( [random.choice(sample) for _ in range(len(sample))] for _ in range(n) ) ) )
        return [ np.quantile( replica_values, q ) for q in quantiles ]

if __name__ == "__main__":

    n = 1000
    sample = [ random.gauss(0,1) for i in range(1000) ]
    mean   = lambda data:sum(data)/float(len(data))

    quantiles = boot_strapped_quantiles( sample, mean, 1000 )
      
