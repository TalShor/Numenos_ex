import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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
    plt.savefig("scatter_plot.png")
    plt.show()
    plt.figure(figsize=(8,6))