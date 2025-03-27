import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from stats import run_all_reports_bc

def antibody_cytoscore_plot(df: pd.DataFrame) -> None:
    """
    Plot a scatter plot of the count of antibodies reads vs cyto score.
    The points are colored based on the status (Pre or On) and arrows indicate the change from Pre to On for each patient.

    Parameters:
    df (pd.DataFrame): DataFrame containing the data to plot. It should have columns 'count_antibodies_reads', 'cyto_score', 'status', and 'pt'.
    """

    corr = df[['count_antibodies_reads', 'cyto_score']]\
        .assign(log_count_antibodies_reads = lambda df: np.log(df.count_antibodies_reads))\
        .corr().loc['log_count_antibodies_reads', 'cyto_score']

    plt.figure(figsize=(8,6))
    # Scatter plot: red dots for Pre, blue dots for On.
    pre = df[df['status'] == 'Pre']
    on = df[df['status'] == 'On']
    plt.scatter(pre['cyto_score'], pre['count_antibodies_reads'], color='red', label='Pre', s=50)
    plt.scatter(on['cyto_score'], on['count_antibodies_reads'], color='blue', label='On', s=50)

    # Draw arrows from Pre to On for matching pt.
    for pt in df['pt'].unique():
        pre_row = df[(df['pt'] == pt) & (df['status'] == 'Pre')]
        on_row = df[(df['pt'] == pt) & (df['status'] == 'On')]
        if not pre_row.empty and not on_row.empty:
            x_start = pre_row['cyto_score'].iloc[0]
            y_start = pre_row['count_antibodies_reads'].iloc[0]
            x_end = on_row['cyto_score'].iloc[0]
            y_end = on_row['count_antibodies_reads'].iloc[0]
            plt.annotate('', xy=(x_end, y_end), xytext=(x_start, y_start),
                    arrowprops=dict(arrowstyle="->", color='gray', lw=3))
            # Add sample name annotation to Pre nodes
            plt.annotate(pre_row['pt'][0], (x_start, y_start), 
                        textcoords="offset points", xytext=(5,5), ha='left')

    plt.yscale('log')
    plt.xlabel("Cyto Score")
    plt.ylabel("Count of Post-filtering Antibodies Reads")
    plt.legend()
    plt.title("Antibodies Reads vs Cyto Score\nCorrelation: {:.2f} of log(#reads) to Cyto Score".format(corr))
    plt.tight_layout()
    plt.show()
    


def bc_heatmap(report_paths:dict) -> None:
    """
    Plot a heatmap of the normalized Bray-Curtis dissimilarity matrix.
    The diagonal of the heatmap is highlighted in red.

    Parameters:
    report_paths (dict): A dictionary containing the paths to the trust report files.
    """

    all_bc = run_all_reports_bc(report_paths)
    normed_all_bc = (all_bc - all_bc.min().min()) / (all_bc.max().max() - all_bc.min().min())

    sns.heatmap(normed_all_bc, annot=True, fmt=".2f", cmap="YlGnBu", cbar=False)

    plt.title("Heatmap of Normalized Bray-Curtis Dissimilarity")
    for i in range(0, normed_all_bc.shape[0], 2):
        plt.gca().add_patch(plt.Rectangle((i, i), 2, 2, linewidth=2, edgecolor='red', facecolor='none'))
        

def plot_max_interactions_boxplot(pairs_dist: pd.DataFrame, n_top_genes=30) -> None:
    """
    Plot a boxplot of the top n VJ gene pairs based on their frequency.
    The boxplot shows the distribution of frequencies for each VJ gene pair.
    
    Parameters:
        pairs_dist (pd.DataFrame): DataFrame containing the frequency of V and J gene pairs for each sample.
        n_top_genes (int): The number of top VJ gene pairs to plot. Default is 30.
    """
    top_pairs_dist = pairs_dist[pairs_dist.median().sort_values(ascending=False).head(n_top_genes).index]
    top_pairs_dist = top_pairs_dist.melt(ignore_index=False).assign(VJ_gene = lambda x: x['V_gene'] + '*' + x['J_gene'])

    plt.figure(figsize=(12, 8))
    sns.boxplot(top_pairs_dist, x='VJ_gene', y='value')
    _ = plt.xticks(rotation=90)


def plot_max_interactions_clustermap(pairs_dist: pd.DataFrame, n_top_genes=100) -> None:
    """
    Plot a heatmap of the top n VJ gene pairs based on their frequency.
    The heatmap shows the frequency of each VJ gene pair for each sample.
    
    Parameters:
        pairs_dist (pd.DataFrame): DataFrame containing the frequency of V and J gene pairs for each sample.
        n_top_genes (int): The number of top VJ gene pairs to plot. Default is 30.
    """
    r = pairs_dist.median().to_frame().sort_values(by=0, ascending=False).head(n_top_genes)\
        .pivot_table(index='V_gene', columns='J_gene', values=0,).fillna(0)

    sns.clustermap(r, xticklabels=True, yticklabels=True, cmap="YlGnBu", cbar=True)