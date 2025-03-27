# TRUST4 Output Files Guide

## Overview

TRUST4 (TCR and BCR Repertoire Utilities for Solid Tissue) is a computational tool that reconstructs immune receptor repertoires from RNA-seq data. When TRUST4 is run, it produces several output files that contain information about the assembled immune receptor sequences. This document explains the purpose and structure of each output file.

## Main Output Files

### 1. Final and Raw Assemblies

**Files:** `TRUST_example_final.out`, `TRUST_example_raw.out`

These files contain the assembled contigs and their coverage information:

```
>assemble8 IGLV2-23+IGLV2-11+IGLV2-14
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3...
```

- **Header format:** `>assembleX [gene names]`
- **Data rows:** Numbers representing read coverage at each position
- `raw.out` contains initial assemblies
- `final.out` contains refined assemblies after additional processing

### 2. Annotated Assemblies

**File:** `TRUST_example_annot.fa`

FASTA format file containing annotated consensus assemblies:

```
>assemble44 333 4.13 TRAV12-2(285):(0-285):(0-285):100.00,TRAV12-1(285):(0-285):(0-285):94.39 +(0):0.00 TRAJ29(31):(285-316):(0-31):100.00 TRAC(295):(316-611):(0-295):100.00 CDR1(64-84):100.00=ATATTACAGGAGGGGGTAAAGC CDR2(132-141):100.00=AAGACTCCG CDR3(277-307):100.00=TGTGCTACGGACCCAACCAGGGATGGATAGCAGCTAT
GAAGAGGATGTTAACATGACTTTGCTTGGATCTGGAACTCCGGTTTCTATTCTGTGTGCTATGAGGTGGGGAAAGCAGGGAGCATCCTGAACTGCACTTATGAGAACAGTGCTTTTGACTACTTCATTTGGTACAGGCAGCGCCTGGGGCAGGGCC...
```

- Each header includes:
  - Consensus ID (`assembleX`)
  - Sequence length
  - Average coverage
  - Gene annotations (V, D, J, C regions)
  - CDR1, CDR2, and CDR3 information with coordinates and sequences

### 3. CDR3 Details 

**File:** `TRUST_example_cdr3.out`

Tab-separated values file with detailed information about the CDR3 regions:

```
consensus_id	index_within_consensus	V_gene	D_gene	J_gene	C_gene	CDR1	CDR2	CDR3	CDR3_score	read_fragment_count	CDR3_germline_similarity	complete_vdj_assembly
assemble440	0	IGHV2-5*01,IGHV2-70*12	IGHD6-13*01	IGHJ4*02	*	*	*	TGTGCACACACTAATGCGGCTATAACAGCAGCAGAATCATTTGACTACTGG	1.00	11.00	96.88	0
```

- **CDR3_score:** 1.00 is maximum, 0.01 means imputed CDR3
- **germline_similarity:** Percentage similarity to germline sequence
- **complete_vdj_assembly:** 1 if complete, 0 if partial

### 4. Simplified Report

**File:** `TRUST_example_report.tsv`

A user-friendly summary focusing on CDR3:

```
read_count	frequency	CDR3_dna	CDR3_amino_acids	V	D	J	C	consensus_id	consensus_id_complete_vdj
11	0.065	TGTGCACACACTAATGCGGCTATAACAGCAGCAGAATCATTTGACTACTGG	CATHNSAITSSSIHLTTY	IGHV2-5	IGHD6-13	IGHJ4	*	assemble440	0
```

- **frequency:** Proportion of reads (BCR/IGH and TCR/TR chains normalized separately)
- **CDR3_amino_acids:** Translated CDR3 sequence, with "_" for stop codons and "?" for ambiguous nucleotides

### 5. AIRR Format Files

**Files:** `TRUST_example_airr.tsv`, `TRUST_example_airr_align.tsv`

Files following the Adaptive Immune Receptor Repertoire (AIRR) standardized format:
- Contains standardized fields for immune repertoire data
- Compatible with other AIRR-compliant analysis tools

## Intermediate Files

### 6. Assembled Reads

**File:** `TRUST_example_assembled_reads.fa`

Contains the actual sequence reads that were assembled into consensus sequences.

### 7. Candidate Reads

**Files:** `TRUST_example_toassemble_1.fq`, `TRUST_example_toassemble_2.fq`

FASTQ files containing reads identified as potential immune receptor sequences:
- These are the reads extracted from the input data that TRUST4 used for assembly
- Used as input for the de novo assembly process

## Interpreting Results

### Nomenclature
- **IGH/IGK/IGL:** Immunoglobulin Heavy/Kappa/Lambda chains (B cells)
- **TRA/TRB/TRG/TRD:** T-cell Receptor Alpha/Beta/Gamma/Delta chains
- **V/D/J/C:** Variable/Diversity/Joining/Constant gene segments

### CDR3 Information
The CDR3 (Complementarity-Determining Region 3) is the most variable part of the receptor and crucial for antigen recognition. TRUST4 focuses on accurate assembly and annotation of CDR3 regions as they are the most important determinants of receptor specificity.

### Coverage Analysis
The numeric values in the raw and final output files represent read coverage at each position, providing confidence metrics for the assembled sequences.