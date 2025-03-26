import pandas as pd
import numpy as np
from scipy.spatial.distance import braycurtis
from read_files import read_trust_report_file_with_header

def report_braycurtis(rep1:pd.DataFrame, rep2:pd.DataFrame) -> float:
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


def run_all_reports_bc(reports:dict) -> pd.DataFrame:
    """
    Calculate the Bray-Curtis distance for all pairs of repertoires in the reports dictionary.
    
    Parameters:
        reports (dict): A dictionary containing the paths to the trust report files.
        The keys are the sample names, and the values are the file paths.

    Returns:
        pd.DataFrame: A DataFrame containing the Bray-Curtis distances between all pairs of repertoires.
    """

    all_reports = {key:read_trust_report_file_with_header(path) for key, path in reports.items()}

    all_bc = []
    for key1, value1 in all_reports.items():
        for key2, value2 in all_reports.items():
            if key1 != key2:
                all_bc.append((key1, key2, report_braycurtis(value1, value2)))
    
    all_bc = pd.DataFrame(all_bc, columns=["Sample1", "Sample2", "Bray-Curtis"])\
        .pivot_table(index="Sample1", columns="Sample2", values="Bray-Curtis")

    return all_bc
    