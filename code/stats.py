import pandas as pd
import numpy as np
from scipy.spatial.distance import braycurtis
from read_files import read_trust_report_file_with_header, get_cyto_score, get_cyto_score_extended
from typing import List, Dict
from scipy.stats import wilcoxon


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
    

def get_file_stats(srr_prefix: str) -> List[int]:
    """
    Get the number of lines in the file with the given prefix.

    Parameters:
        srr_prefix (str): The prefix of the file name.

    Returns:
        List[int]: A list containing the number of lines in the file.
    """

    def _get_file_length(file_path: str) -> pd.Series:
        """
        Get the number of lines in a file.

        Parameters:
            file_path (str): The path to the file.

        Returns:
            int: The number of lines in the file.
        """
        with open(file_path, "r") as f:
            return sum(1 for _ in f)
    
    results = {}
    results['antibody_raw_count'] = _get_file_length(srr_prefix + "_raw.out") / 6
    results['antibody_final_count'] = _get_file_length(srr_prefix + "_final.out") / 6
    results['total_read_count'] = _get_file_length(srr_prefix + "_1.fastq") / 4
    return pd.Series(results)


def get_file_stats_from_list(report_paths: dict, antibody_reads:pd.DataFrame) -> pd.DataFrame:
    results = {
        name:get_file_stats(path.replace('_report.tsv', '')) \
            for name, path in report_paths.items()
    }
    return pd.DataFrame(results).T.join(antibody_reads)


def get_vdjc_aggregated(reports_dfs: dict) -> pd.DataFrame:
    """
    Aggergate the V(D)J(C) genes from the trust report files.

    Parameters:
        report_paths (dict): A dictionary containing the paths to the trust report files.
        The keys are the sample names, and the values are the file paths.
    Returns:
        pd.DataFrame: A DataFrame containing the V(D)J(C) genes for each sample.
    """
    return pd.concat(reports_dfs.values).set_index(['sample', 'frequency'])\
        [['V', 'D', 'J', 'C']].melt(ignore_index=False).reset_index()\
        .groupby(['sample', 'variable', 'value'])['frequency'].sum().sort_values()\
        .reset_index().pivot_table(
            index=['sample'], 
            columns=['variable', 'value'], 
            values='frequency', fill_value=0)

def get_vj_genes_freq(reports_dfs: dict) -> pd.DataFrame:
    """
    Aggregate the V and J genes from the trust report files.
    Parameters:
        reports_dfs (dict): A dictionary containing the DataFrames of the trust report files.
        The keys are the sample names, and the values are the DataFrames.
    Returns:
        pd.DataFrame: A DataFrame containing the frequency of V and J genes for each sample.
    """

    return pd.concat(reports_dfs.values)\
        .assign(V_gene = lambda x: x['V'].str.split('*').str[0])\
        .assign(J_gene = lambda x: x['J'].str.split('*').str[0])\
        .groupby(['sample', 'V_gene', 'J_gene'])[['frequency']].sum()\
        .sort_values(by='frequency', ascending=False)

def get_vj_genes_pairs_dist(VJ_genes: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate the distribution of V and J gene pairs.
    Parameters:
        VJ_genes (pd.DataFrame): A DataFrame containing the frequency of V and J genes for each sample.
    Returns:    
        pd.DataFrame: A DataFrame containing the distribution of V and J gene pairs.
    """
    return VJ_genes.reset_index().query('V_gene != "." and J_gene != "."')\
        .pivot_table(
            index='sample', 
            columns=['V_gene', 'J_gene'], 
            values='frequency', 
            fill_value=0)

def compare_genes_thought_test(genes:pd.DataFrame, gene_col_name:str) -> pd.DataFrame:
    """
    Compare the frequencies of genes in the 'On' and 'Pre' conditions using the Wilcoxon test.
    Parameters:
        genes (pd.DataFrame): A DataFrame containing the frequencies of genes.
        gene_col_name (str): The name of the column containing the gene names.
    Returns:
        pd.DataFrame: A DataFrame containing the p-values of the Wilcoxon test for each gene.    
    """

    def run_wilcoxon(df):
        df = df.T
        try:
            on_df = df.loc[:, (slice(None), "On")]
            pre_df = df.loc[:, (slice(None), "Pre")]
            return wilcoxon(on_df, pre_df)
        except:
            return None

    cyto_score = get_cyto_score_extended(get_cyto_score())
    genes_time_compare = genes.join(cyto_score.reset_index(level=1)[['pt', 'status']])\
        .pivot_table(index='pt', columns=[gene_col_name, 'status'], values='frequency', aggfunc='sum').fillna(0)
        
    genes_time_compare_pvalue = genes_time_compare.T.groupby(gene_col_name)\
        .apply(run_wilcoxon).dropna()
    
    genes_time_compare_pvalue = genes_time_compare_pvalue.to_frame('statistic')\
        .assign(pvalue=lambda df: df.statistic.apply(lambda x: x.pvalue[0])).sort_values(by='pvalue').drop(columns='statistic')
    return genes_time_compare_pvalue