from pathlib import Path
import pandas as pd

def read_cdr3_out_file(path:Path) -> pd.DataFrame:
    """
    Read the CDR3 output file from the specified path.
    The file is expected to be a tab-separated values (TSV) file with the following columns:
    - consensus_id
    - index_within_consensus
    - V_gene
    - D_gene
    - J_gene
    - C_gene
    - CDR1
    - CDR2
    - CDR3
    - CDR3_score
    - read_fragment_count
    - CDR3_germline_similarity
    - complete_vdj_assembly

    Parameters:
        path (Path): The path to the CDR3 output file.
    Returns:
        pd.DataFrame: A DataFrame containing the CDR3 output data.
    """

    return pd.read_csv(path, sep='\t', 
        names=['consensus_id', 'index_within_consensus', 
               'V_gene', 'D_gene', 'J_gene', 'C_gene', 
               'CDR1', 'CDR2', 'CDR3', 'CDR3_score', 
               'read_fragment_count', 'CDR3_germline_similarity', 
               'complete_vdj_assembly'], header=None)



def read_trust_report_file(path:Path) -> pd.DataFrame:
    """
    Read the trust report file from the specified path.
    The file is expected to be a tab-separated values (TSV) file with the following columns:
    - read_count
    - frequency
    - CDR3_dna
    - CDR3_amino_acids
    - V
    - D
    - J
    - C
    - consensus_id
    - consensus_id_complete_vdj

    Parameters:
        path (Path): The path to the trust report file.
    Returns:
        pd.DataFrame: A DataFrame containing the trust report data.
    """

    return pd.read_csv(path, sep='\t', 
        names=['read_count', 'frequency', 'CDR3_dna', 'CDR3_amino_acids', 
               'V', 'D', 'J', 'C', 'consensus_id', 'consensus_id_complete_vdj'], header=None)