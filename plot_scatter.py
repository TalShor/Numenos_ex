import matplotlib.pyplot as plt
import pandas as pd
from io import StringIO

# Example input data; adjust to your file location if needed
data = """count_antibodies_reads pt_sample_name cyto_score pt status
19702 Pt1_Pre_AD101148-6 48.985749 Pt1 Pre
1933 Pt1_On_AD174047-6 13.230992 Pt1 On
1424 Pt2_Pre_AD101150-6 40.385277 Pt2 Pre
529 Pt2_On_AD174046-6 44.807463 Pt2 On
"""

# Read data into a DataFrame (modify separator as needed for your file)
df = pd.read_csv(StringIO(data), sep=r"\s+")
df['count_antibodies_reads'] = pd.to_numeric(df['count_antibodies_reads'])
df['cyto_score'] = pd.to_numeric(df['cyto_score'])

plt.figure(figsize=(8,6))
# Scatter plot: red dots for Pre, blue dots for On.
pre = df[df['status'] == 'Pre']
on = df[df['status'] == 'On']
plt.scatter(pre['cyto_score'], pre['count_antibodies_reads'], color='red', label='Pre')
plt.scatter(on['cyto_score'], on['count_antibodies_reads'], color='blue', label='On')

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
                     arrowprops=dict(arrowstyle="->", color='gray', lw=1.5))

plt.xlabel("cyto_score")
plt.ylabel("count_antibodies_reads")
plt.legend()
plt.title("Scatter Plot: Antibodies Reads vs cyto_score")
plt.tight_layout()
plt.savefig("scatter_plot.png")
plt.show()
