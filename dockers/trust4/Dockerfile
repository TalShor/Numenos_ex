FROM continuumio/miniconda3

RUN conda install conda-forge::mamba

RUN mamba install bioconda::trust4

COPY hg38_bcrtcr.fa /reference/hg38_bcrtcr.fa
COPY human_IMGT+C.fa /reference/human_IMGT+C.fa