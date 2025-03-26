import pandas as pd
import numpy as np
from scipy.spatial.distance import braycurtis

def braycurtis(rep1:pd.DataFrame, rep2:pd.DataFrame) -> float:
    """
    Calculate the Bray-Curtis distance between two repertoires.
    
    Parameters:
        rep1 (pd.DataFrame): The first repertoire DataFrame.
        rep2 (pd.DataFrame): The second repertoire DataFrame.

    Returns:
        float: The Bray-Curtis distance between the two repertoires.
    """
    
    compare_freqs = rep1.set_index('CDR3aa')[['frequency']]\
        .join(rep2.set_index('CDR3aa')[['frequency']], rsuffix='_rep2', lsuffix='_rep1', how='outer').fillna(0)\
        .drop('out_of_frame')
    
    return braycurtis(compare_freqs['frequency_rep1'], compare_freqs['frequency_rep2'])
