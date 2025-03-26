main project 
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061

the name of the paper is 
Molecular portraits of tumor mutational and micro-environmental sculpting by immune checkpoint blockade therapy

main ftp site - 
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/

cyto score
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061_BMS038109Sample_Cytolytic_Score_20161026.txt.gz

Key Immunophenotyping File
GSE91061_BMS038109Sample_Cytolytic_Score_20161026.txt.gz

This contains the Cytolytic Score for each sample
The cytolytic score is a direct measure of immune cell killing activity in tumors
This score is calculated from the expression of cytolytic genes (GZMA and PRF1)
Higher scores = "hot" tumors with active immune infiltration
Lower scores = "cold" tumors with limited immune activity
This file is perfect for selecting samples based on immune status. Here's how to use it:

translate_GSM_to_pt.csv 
translate between GSM name, pt name for the experiment and the SRR name 
we download the data via prefetch from the srr-tools