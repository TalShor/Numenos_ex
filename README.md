# Numenos_ex Project

## Project of our dataset
[Molecular portraits of tumor mutational and micro-environmental sculpting by immune checkpoint blockade therapy](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061)

### Paper 
[Tumor and Microenvironment Evolution during Immunotherapy with Nivolumab](https://pmc.ncbi.nlm.nih.gov/articles/PMC5685550/)

The paper investigates the genomic and microenvironmental changes in melanoma tumors during anti-PD-1 immunotherapy with nivolumab. It aims to understand how immune checkpoint blockade influences tumor mutation load, clonal evolution, and *immune cell dynamics*. By analyzing genomic, *transcriptomic*, and *T-cell receptor data*, the study provides insights into the mechanisms of response and resistance to immunotherapy.

### Data
[Downloadalbe via SRA tools](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA356761)

Contains 
* 68 pairs of RNA-Seq (before treatment with Nivolumab and after) for PBMCs (50M reads post QC)
* Cyto Score (GSE91061_BMS038109Sample_Cytolytic_Score_20161026.txt) based on RAN expression levels

## Preprocessing

### Tools
Orchestration done with [snakemake](./ETL/download_and_preprocess.smk) and job-specific [dockers](./dockers/)

#### ETL Flow
1. Download SRA from NCBI
2. Break SRA apart into .fastq
3. Run quality assesment (FastQC)
4. Optional: Run fastp for trimming (commented out this part because the files are already post QC and in a high quality - see [SRR5088813](./analysis/SRR5088813_1_fastqc.html) for example)
5. Run TRUST4 over the 2 fastq. Use hg38 and TRUST4's references

## Analysis

### [Sample Analysis](./analysis/Sample_Analysis.ipynb)

Testing the reportoires of individuals compared to their counterpart test (Pre / On Nivolumab) and their peers. (using the CDR3aa for this part)

#### Outcomes
1. A reportoire of an individual is fairy consistent even if taking Nivolumab (at least it's more similar than other individuals samples) 
2. Apart from SRR5088829 - probably due to a significant low number of V(D)J aligned reads


### [Clonotype Diversity](./analysis/Clonotype_Diversity.ipynb)

After aggregating differnet variants for same genes for the V/J parts 
Showing interactions between V/J genes
Showing top visible genes for each of the families

#### Outcomes

1. Some interactions are common between most samples
2. Huge flactuation in percentages between samples for the same genes

