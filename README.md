# Numenos.ai Home Task

## Table of Contents
1. [Project of our dataset](#project-of-our-dataset)
2. [Paper](#paper)
3. [Data](#data)
4. [Preprocessing](#preprocessing)
    - [Tools](#tools)
    - [ETL Flow](#etl-flow)
    - [Usage](#usage)
5. [Analysis](#analysis)
    - [Sample Analysis](#sample-analysis)
    - [Clonotype Diversity](#clonotype-diversity)
6. [Architecture](#architecture)
7. [Challenges and Solutions](#challenges-and-solutions)
    - [Technical](#technical)
    - [Scientific](#scientific)
8. [Reliability of the output discussion](#reliability-of-the-output-discussion)

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

### ETL Flow
1. Download SRA from NCBI
2. Break SRA apart into .fastq
3. Run quality assesment (FastQC)
4. Optional: Run fastp for trimming (commented out this part because the files are already post QC and in a high quality - see [SRR5088813](./analysis/SRR5088813_1_fastqc.html) for example)
5. Run TRUST4 over the 2 fastq. Use hg38 and TRUST4's references

### Usage
`snakemake -s ETL/download_and_preprocess.smk`

### Requirements
1. docker installed
2. run all the build scripts in the dockers directry

## Analysis

### [Sample Analysis](./analysis/Sample_Analysis.ipynb)

Testing the reportoires of individuals compared to their counterpart test (Pre / On Nivolumab) and their peers. (using the CDR3aa for this part)

#### Outcomes
1. A reportoire of an individual is fairy consistent even if taking Nivolumab (at least it's more similar than other individuals samples) 
2. Apart from SRR5088829 - probably due to a significant low number of V(D)J aligned reads (and maybe SRR5088830? unsure)

### [Clonotype Diversity](./analysis/Clonotype_Diversity.ipynb)

After aggregating differnet variants for same genes for the V/J parts 
Showing interactions between V/J genes
Showing top visible genes for each of the families

#### Outcomes

1. Some interactions are common between most samples
2. Huge flactuation in percentages between samples for the same genes
3. The V / J pair IGKV1-39 / IGKJ3 might have some significance for differentiating between Pre-medication and On-medication

### Requirements
numpy, pandas, scipy, matplotlib, seaborn, 

## Architecture

1. [analysis](./analysis/) - Jupyetr notebooks for the analysis
2. [code](./code/) - everything from plotting, transformation of data and statistics
3. [data](./data/) - not in .github. where we store all the pre-and post processed data
4. [dockers](./dockers) - dockers used in the ETL
5. [documents](./documents/) - papers, technical documentations etc.
6. [ETL](./ETL/) - snakemake code 


## Challenges and Solutions

### Technical
1. Getting the Data  
TCGA is kept tightly under wraps, so it was surprising to find a study with so much freely downloadable information. It took quite a while to identify one with good phenotypes accompanied by PBMC RNA-Seq data.  
**Solution:** I used the GEO access tool suggested in the exercise instead of TCGA.

2. Working with Snakemake on macOS  
It took longer than expected due to some unsupported flags. I had worked with Snakemake a couple of years ago on Linux, so some adjustments were necessary.  
**Solution:** Following guidance from ChatGPT, I reconfigured a couple of settings to accommodate the differences between Linux and macOS for Snakemake and Docker.

### Scientific 
1. New Field  
I had experience with the protein side of immunoassays but not with antibodies themselves. The jargon was new, and the data required some formatting (e.g., aggregating all variants of the genes).  
**Solution:** I read the TRUST paper, mapped the output according to the [README](https://github.com/liulab-dfci/TRUST4/blob/master/README.md) and experimented with the data before composing the analysis.

2. Sample Similarity  
Spearman and Pearson correlations did not provide meaningful results due to a zero-inflated distribution when comparing two samples. Jaccard removed much of the information by being binary.  
**Solution:** I ended up using Bray-Curtis distance. Since the range was very small (with similar samples around ~0.97 BC distance and the maximum at 1.0), I normalized the distances to enhance visibility in the graph.

3. Hard to show significance
* The gene frequencies are zero tailed (the mapping is sparse)
* It's log-normal at nature (probably due to some cascading effect) when not missing alltogether.
* *The numbers are very small
We can't use a lot of tests (e.g. even though chi-square works for small numbers - it needs a normal distribution)
**Solution:** Use Wilcoxon test. Doesn't assume normal distribution and can operate with a small number of samples.
That being said - the results were lackluster with 7 samples (Pre and On) and if we would have corrected for multiple testing - we would have been left with nothing.


## Reliability of the output discussion

1. Run FastQC to validate that the reads are of a sufficient quality.
2. Use the sum of the frequencies to easily detect if all types (families) are accounted for
3. Check nubmer of assembled sequences (in the final.out)
4. Optional: Since it's RNA-Seq of PBMCs - we can check B-cells / T-cells ratio as it should be somewhat correlative. 