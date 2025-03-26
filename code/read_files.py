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


def read_trust_report_file_with_header(path:Path) -> pd.DataFrame:
    """
    Read the trust report file from the specified path.
    The file is expected to be a tab-separated values (TSV) file with the following columns:
    - #count
    - frequency
    - CDR3nt
    - CDR3aa
    - V
    - J
    - C
    - cid
    - cid_full_length

    Parameters:
        path (Path): The path to the trust report file.
    Returns:
        pd.DataFrame: A DataFrame containing the trust report data.
    """

    return pd.read_csv(path, sep='\t')


def get_cyto_score() -> pd.DataFrame:
    """
    Read the cytolytic score data from the specified files.
    The function reads two files:
    - 'data/translate_GSM_to_pt.csv': A CSV file containing sample names.
    - 'GSE91061_BMS038109Sample_Cytolytic_Score_20161026.txt': A tab-separated file containing cytolytic scores.
    The function joins the two DataFrames on the sample names and returns the resulting DataFrame.
    Returns:
        pd.DataFrame: A DataFrame containing the cytolytic scores.
    """
    sample_names = pd.read_csv('../data/translate_GSM_to_pt.csv')\
        .set_index('pt_sample_name')
    cyto_scores = pd.read_csv('../data/GSE91061_BMS038109Sample_Cytolytic_Score_20161026.txt', sep='\t')\
        .set_index('V1').rename(columns={'V2':'cyto_score'})

    return sample_names.join(cyto_scores, how='inner')\
        .dropna(subset=['srr_sample_name'])\
        .reset_index().rename(columns={'index':'pt_sample_name'})\
        .set_index(['srr_sample_name', 'pt_sample_name'])\
        .drop(columns=['gsm_sample_name'])


def get_paried_cyto_score(cyto_score:pd.DataFrame) -> pd.DataFrame:
    info = cyto_score.index.get_level_values('pt_sample_name').str.split('_')
    cyto_score['pt'] = info.str[0]
    cyto_score['status'] = info.str[1]
    cyto_score_pv = cyto_score.pivot_table(
        index=['pt'], columns=['status'], values='cyto_score')
    
    return cyto_score_pv[['Pre', 'On']]