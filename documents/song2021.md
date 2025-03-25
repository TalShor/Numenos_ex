<img src="./fo4fqvbb.png"
style="width:0.98264in;height:0.18576in" /><img src="./lbwelckf.png" style="width:0.7092in" /><img src="./kt4nwace.png"
style="width:0.13368in;height:0.13281in" /><img src="./jkjz4oo2.png"
style="width:1.39811in;height:0.2934in" /><img src="./prdts0lp.png"
style="width:1.07118in;height:0.23872in" />

**TRUST4:** **immune** **repertoire** **reconstruction** **from**
**bulk** **and** **single-cell** **RNA-seq** **data**

**Li** **Song1,2,** **David** **Cohen1,** **Zhangyi** **Ouyang **
**[ ](http://orcid.org/0000-0002-9905-9832)3,** **Yang** **Cao4,**
**Xihao** **Hu1,2** **and** **X.** **Shirley**
**Liu**[** **](http://orcid.org/0000-0003-4736-7339) ** 1,2,5** ✉

**We** **introduce** **the** **TRUST4** **open-source** **algorithm**
**for** **recon-struction** **of** **immune** **receptor**
**repertoires** **in** **αβ/γδ** **T cells** **and** **B cells**
**from** **RNA-sequencing** **(RNA-seq)** **data.** **Compared**
**with** **competing** **methods,** **TRUST4** **supports** **both**
**FASTQ** **and** **BAM** **format** **and** **is** **faster** **and**
**more** **sensitive** **in** **assembling** **longer—** **even**
**full-length—receptor** **repertoires.** **TRUST4** **can** **also**
**call** **repertoire** **sequences** **from** **single-cell**
**RNA-seq** **(scRNA-seq)** **data** **without** **V(D)J**
**enrichment,** **and** **is** **compatible** **with** **both**
**SMART-seq** **and** **5′ 10x Genomics** **platforms.**

Both T and B cells can generate diverse receptor (TCR and BCR,
respectively) repertoires, through somatic V(D)J recombi-nation, to
recognize various external antigens or tumor neoanti-gens. Following
antigen recognition, BCRs also undergo somatic hypermutations (SHMs) to
further improve antigen-binding affinity. Repertoire sequencing has been
increasingly adopted in infectious disease1, allergy2, autoimmune3,
tumor immuology4 and cancer immunotherapy5 studies, but it is an
expensive assay and consumes valuable tissue samples. Alternatively,
RNA-seq data contain expressed TCR and BCR sequences in tissues or
peri-pheral blood mononuclear cells (PBMC). However, because reper-toire
sequences from V(D)J recombination and SHM are different from the
germline, they are often eliminated in the read-mapping step.

Previously we developed the TRUST algorithm6–8, utilized to de novo
assemble immune receptor repertories directly from tis-sue or blood
RNA-seq data. When applied to The Cancer Genome Atlas (TCGA) tumor
RNA-seq data, TRUST revealed profound biological insights into the
repertoires of tumor-infiltrating Tcells6 and Bcells8, as well as their
associated tumor immunity. Although less sensitive than TCR-seq and
BCR-seq, TRUST is able to iden-tify the abundantly expressed and
potentially more clonally expanded TCRs/BCRs in the RNA-seq data that
are more likely to be involved in antigen binding9. Recent years have
also seen other computational methods introduced for immune repertoire
construction from RNA-seq data, including V’DJer10, MiXCR11, CATT12 and
ImRep13. These methods focus on reconstruction of
complementary-determining region3 (CDR3), with limited ability to
assemble full-length V(D)J receptor sequences, although CDR1 and CDR2 on
the Vsequence still contribute considerably to anti-gen recognition and
binding. For example, five out of six mutations predicted in a recent
study to influence antibody affinity in the acidic tumor environment are
located in CDR1 and CDR2 (ref. 14), and four out of nine positions
contributing most to 4A8 antibody binding to the SARS-CoV-2 spike
protein are in CDR1 and CDR2 (ref. 15). Therefore, algorithms that can
infer full-length immune receptor repertoires can facilitate better
receptor–antigen interac-tion modeling.

With the advance of scRNA-seq technologies, researchers can study immune
cell gene expression and receptor repertoire sequences simultaneously.
Several algorithms, including MiXCR11, BALDR16, BASIC17 and VDJPuzzle18,
have been developed to con-struct full-length paired TCRs or BCRs from
the SMART-seq scRNA-seq platform19. In contrast to SMART-seq,
droplet-based scRNA-seq platforms such as 10xGenomics20, while yielding
sparser transcript coverage per cell, can process orders of magni-tude
more cells at lower cost. To analyze immune repertoires using the
10xGenomics platform, researchers currently need to prepare extra
libraries to amplify TCR/BCR sequences.

In this study, we redesigned the TRUST algorithm to TRUST4 with
substantially enhanced features and improved performance for immune
repertoire reconstruction (Fig. 1a). First, TRUST4 supports fast
extraction of TCR/BCR candidate reads from either FASTQ or BAM files.
Second, TRUST4 prioritizes candidate read assembly by abundance and
assembles all candidate reads with par-tial overlaps against contigs,
thus increasing algorithm speed. Third, TRUST4 explicitly represents
highly similar reads in the contig con-sensus, thus accommodating
somatic hypermutations and improv-ing memory efficiency (Methods).
Fourth, TRUST4 can assemble full-length V(D)J sequences on TCRs and
BCRs. Finally, TRUST4 supports repertoire reconstruction from scRNA-seq
platforms without requiring the extra 10x V(D)J amplification steps.

We evaluated the performance of TRUST4 on TCR/BCR recon-struction from
bulk RNA-seq using three different approaches. First, for TCR evaluation
we used in silico RNA-seq datasets with known TRB sequences from an
earlier study11. On average, TRUST4 called 281% more CDR3s than MiXCR,
22.9% more than CATT, 57.8% more than TRUST3 and maintained a zero
false-positive rate across different read lengths (Fig. 1b; further
parameter set-tings are given in Supplementary Fig. 1a). Second, for BCR
evalu-ation, we used six tumor RNA-seq samples of ~100million pairs of
150base pair (bp) reads with corresponding immunoglobulin heavy-chain
(IGH) BCR-seq as the gold standard8. Since BCRs also have somatic
hypermutation and isotype switching during clonal expansion, we required
the algorithm call to match CDR3 and V, J and C (isotype) gene
assignments as BCR-seq. TRUST4 showed better precision (\>18%) and
sensitivity (\>74%) than MiXCR in five out of six samples (Fig. 1c;
further parameter settings are given in Supplementary Fig. 2a). On the
sixth sample, TRUST4 lost only 6% precision with twice the sensitivity
compared to MiXCR (FZ-97). We note that BCR-seq and RNA-seq were
conducted on dif-ferent slices of the same tumor. Even two technical
replicates of repertoire sequencing on the same DNA/RNA could not
achieve 100% precision and sensitivity, so the performance metrics are
likely to be underestimations. TRUST4 consistently assembled

1Department of Data Science, Dana-Farber Cancer Institute, Boston, MA,
USA. 2Harvard T.H. Chan School of Public Health, Boston, MA, USA.
3Department of Biotechnology, Beijing Institue of Radiation Medicine,
Beijing, China. 4College of Life Sciences, Sichuan University, Chengdu,
China. 5Center for Functional Cancer Epigenetics, Dana-Farber Cancer
Institute, Boston, MA, USA. e-mail: <xsliu@ds.dfci.harvard.edu>

**NATURe** **MeTHODS**\| VOL 18 \| JUnE2021 \| 627–630 \|
[www.nature.com/naturemethods](http://www.nature.com/naturemethods)
**627**

Brief CommuniCation **Nature** **Methods**

> **a**
>
> <img src="./fsk4kqjd.png"
> style="width:0.73667in;height:0.81833in" /><img src="./ywadri1f.png" style="width:3.795in;height:0.395in" /><img src="./1alpuyog.png"
> style="width:0.99436in;height:0.25694in" />BAM Candidate reads Contigs
>
> Annotations

UTR1 IGHV1-2 CDR3 IGHJ1 IGHG1

> or UTR2 TRBV2 CDR3 TRBJ2-3 TRBC
>
> UTR3 TRGV4 CDR3 TRGJ2 TRGC FASTQ
>
> **b** Paired-end
>
> 350
>
> 300
>
> 250
>
> 200
>
> 150
>
> 100
>
> 50
>
> 0

**c** **d**

> <img src="./iexukds4.png"
> style="width:2.01539in;height:0.42068in" /><img src="./0ywr4rea.png"
> style="width:1.09925in;height:1.57929in" /><img src="./504gb5ap.png"
> style="width:0.61633in;height:1.57929in" /><img src="./jijep2u3.png"
> style="width:0.61633in;height:1.57929in" />50 V CDR3 J TRUST4 MIXCR
> TRUST4 MIXCR TRUST4 MIXCR
>
> <img src="./bpn4yntf.png"
> style="width:1.56424in;height:0.86458in" />45 FZ-83 FZ-20 FZ-116
>
> FZ-122
>
> 40 FZ-94
>
> FZ-97
>
> 35 1620PV
>
> 30
>
> VH_CD19pos 25 TRUST4
>
> MiXCR
>
> 20
>
> 0 2,000 4,000 6,000 8,000
>
> No. of recalls
>
> **Fig.** **1** **\|** **The** **performance** **of** **TRUST4** **on**
> **bulk** **RNA-seq** **data.** **a**–**d**, TRUST4 applied to bulk
> RnA-seq data. **a**, Overview of TRUST4. **b**, number of TRB CDR3s
> reported by MiXCR, CATT, TRUST3 and TRUST4 from in silico RnA-seq
> data. **c**, Precision–recall of six bulk RnA-seq samples using
> BCR-seq results as the gold standard. **d**, Evaluation of full-length
> V, CDR3 and J sequences assembled by TRUST4 and MiXCR on pseudobulk
> RnA-seq by grouping SMART-seq data. Each row represents whether the
> cell’s sequences were recovered in the pseudobulk data.
>
> more IGHs across different abundance ranges reported in BCR-seq
> (Supplementary Fig. 2b), and found twice as many IGHs with a single
> RNA copy than MiXCR. In addition, TRUST4 required only 20–25% of the
> time, on average, needed by MiXCR to process these samples
> (Supplementary Table 1), at \<6GB of memory usage on an eight-thread
> processor. Furthermore, TRUST4 run directly on FASTQ files was notably
> faster than read mapping used to gener-ate BAM files, followed by
> TRUST4 run on BAM files. Third, for base-level, full-length assembly
> evaluation, we created pseudobulk RNA-seq data by randomly selecting
> 25 million read pairs from 137SMART-seq Bcells as a test case. To
> establish a gold standard for BCR calls, we used the 128IGH assemblies
> consistently called by BALDR and BASIC at the single-cell level
> (Supplementary Fig. 3a). TRUST4 and MiXCR correctly identified all
> 128CDR3s and TRUST4 reconstructed 93full-length IGH sequences, while
> MiXCR found only 39 (Fig. 1d). TRUST4 was able to call some BCRs with
> only 5,000randomly sampled read pairs in the SMART-seq data-set
> (18read pairs per chain), and showed higher sensitivity than MiXCR
> across all abundance ranges (Supplementary Fig. 3b). The high
> efficiency of TRUST4 allowed us to characterize the immune repertoire
> in tumor samples, and we identified an association of IgA1 B-cell
> clonal expansion with poor prognosis in colon adeno-carcinoma (COAD)
> from TCGA RNA-seq samples (Methods and Supplementary Fig. 4). We note
> that *IGHA1* overexpression is not associated with survival,
> suggesting that immune repertoire analysis provides additional
> insights into tumor immunity.
>
> Next, we evaluated TRUST4 performance on 5′10x Genomics scRNA-seq data
> on PBMC. For this dataset, the two separately pro-cessed T- and B-cell
> 10x V(D)J libraries served as the gold standard. When considering
> single cells that passed the Seurat21 cell-level quality control,
> TRUST4 made 5,091T- and 1,318B-cell calls (Fig. 2a and Supplementary
> Fig. 5a). Among the CDR3s reported by 10x V(D)J, TRUST4 recovered
> 48.1% (6,035/12,558) of TCR CDR3s and 78.0% (1,946/2,494) of BCR
> CDR3s. The higher sensitivity of TRUST4 on BCR is due to the higher
> expression level of BCR in
>
> **628**

Bcells (Supplementary Fig. 5b). For precision, 94.6% of TCR CDR3s and
98.2% of BCR CDR3s from TRUST4 were identical to 10x V(D) J (Fig. 2b).
Although CellRanger_VDJ was designed for 10x V(D) J data, we tested it
on 5′10x scRNA-seq data, which have the same format. TRUST4 found 78%
more TCR CDR3s and 16% more BCR CDR3s in cells that passed quality
control (Supplementary Fig. 5c). In addition, TRUST4 was over ten times
faster and over twice more memory efficient than CellRanger_VDJ.
Furthermore, TRUST4 also reported 83γδT cells, for which 10x V(D)J
currently does not have a kit to profile. In these data, Seurat did not
anno-tate any γδT cells but rather called 71 out of 83TRUST4-annotated
γδT cells as CD8 Tcells. Close examination of gene expression in these
83cells revealed that they had higher δV and δC gene expres-sion but
lower *CD8A* or *CD8B* expression (Supplementary Fig. 5d), supporting
TRUST4’s annotation of these cells as γδT cells.

We further tested TRUST4 on a 10x Genomics non-small cell lung cancer
(NSCLC) dataset. In this case, TRUST4 called 1,241Tcells and 2,478
Bcells (Supplementary Fig. 6). TRUST4 assembled 142IGH CDR3s out of the
144Seurat-annotated plasma Bcells while 10x V(D)J found only 131. For
these plasma Bcells, TRUST4 also reconstructed full-length paired BCRs
for 104cells in which we observed a high correlation for SHM rate
between IGHs and IGK/IGLs (Fig. 2d; Pearson *r*=0.67, *P*=8×10–15),
sug-gesting coordinated SHMs on two chains during B-cell division.
Furthermore, TRUST4 found more somatic hypermutations on IGH than on
IGK/IGL (*P*\<1×10–10, two-sided Wilcoxon signed-rank test), supporting
the more important role of antibody heavy chain in antigen-binding
affinity.

In summary, TRUST4 is an effective method to infer TCR and BCR
repertoires from bulk RNA-seq or scRNA-seq data. TRUST4 not only has
high efficiency, sensitivity and precision in reconstruc-tion of CDR3s,
but can also assemble full-length immune recep-tor sequences from bulk
RNA-seq data. Furthermore, TRUST4 can reconstruct immune receptor
sequences at the single-cell level, including γδT cells, directly from
5′10x Genomics scRNA-seq data

> **NATURe** **MeTHODS**\| VOL 18 \| JUnE2021 \| 627–630 \|
> [www.nature.com/naturemethods](http://www.nature.com/naturemethods)

**Nature** **Methods**

> **a**
>
> <img src="./atr3ssku.png"
> style="width:2.80667in;height:2.73833in" />10
>
> Brief CommuniCation

<img src="./4zpljgjq.png"
style="width:0.1276in;height:0.1441in" /><img src="./l0yhpy05.png"
style="width:0.16059in;height:0.1224in" /><img src="./t3tcr3d2.png"
style="width:0.20833in;height:0.15191in" /><img src="./fufl5jit.png" /><img src="./l1gomauq.png"
style="width:0.17014in;height:0.21615in" /><img src="./hddr0n0k.png"
style="width:0.22743in;height:0.17535in" /><img src="./fn35yypn.png"
style="width:0.15712in;height:0.2092in" /><img src="./rzskrzr4.png"
style="width:0.10503in;height:0.17535in" /><img src="./wdvaddur.png"
style="width:0.41667in;height:0.46441in" /><img src="./a0lmsv3d.png"
style="width:0.3559in;height:0.3125in" /><img src="./umjaxtee.png"
style="width:0.25521in;height:0.27865in" />**c** 100

> Naive CD4 T cells 5
>
> Others
>
> Treg
>
> 0
>
> CD8 T cells
>
> 95

RestDCs

> Monocytes
>
> 90
>
> –5
>
> –10
>
> B cells
>
> αβ T cells

RestNK 85

> Cor = 0.669
>
> Memory B cells
>
> Naive B cells
>
> γδ T cells
>
> 80
>
> –7.5 –5.0 –2.5 0 2.5 5.0 7.5 10.0 12.5 80 85 90 95 100
>
> UMAP 1 IGHV similarity
>
> **b**
>
> 4,000 1,200
>
> 3,500 1,000
>
> 3,000
>
> 800 2,500
>
> 2,000 600
>
> 1,500 400
>
> 1,000
>
> 500 200
>
> 0 0
>
> TRB TRA IGH IGL/IGK
>
> **Fig.** **2** **\|** **The** **performance** **of** **TRUST4** **on**
> **scRNA-seq** **data.** **a**–**c**, Application of TRUST4 to 5′ 10x
> Genomics scRnA-seq data. **a**, Uniform manifold approximation and
> projection (UMAP) of 5′ 10x Genomics PBMC data. **b**, number of CDR3s
> matched with 10x Genomics V(D)J-enriched library from cells annotated
> by Seurat. **c**, Comparison of TRUST4-assembled V gene similarities
> and reference germline V gene sequences from paired, full-length IGH
> and IGK/IGL assemblies of 5′ 10x Genomics nSCLC data. Each dot
> represents one cell. Teg, regulatory T cells; RestnKs, resting natural
> killer cells; RestDCs, resting dendritic cells; cor, Pearson
> correlation.
>
> without specific 10x V(D)J enrichment libraries. Our results sup-port
> the advantage of the 5′10x Genomics scRNA-seq platform, which not only
> provides gene expression information but also enables computational
> calling of immune repertoires. TRUST4 is available open source at
> [https://github.com/liulab-dfci/TRUST4,](https://github.com/liulab-dfci/TRUST4)
> and could be an important method for tumor immunity and immu-notherapy
> studies.
>
> **Online** **content**
>
> Any methods, additional references, Nature Research report-ing
> summaries, source data, extended data, supplementary infor-mation,
> acknowledgements, peer review information; details of author
> contributions and competing interests; and statements of data and code
> availability are available at
> [https://doi.org/10.1038/](https://doi.org/10.1038/s41592-021-01142-2)
> [s41592-021-01142-2](https://doi.org/10.1038/s41592-021-01142-2).
>
> Received: 18 August 2020; Accepted: 2 April 2021; Published online: 13
> May 2021
>
> **NATURe** **MeTHODS**\| VOL 18 \| JUnE2021 \| 627–630 \|
> [www.nature.com/naturemethods](http://www.nature.com/naturemethods)

**References**

> 1\. Lee, J. et al. Molecular-level analysis of the serum antibody
> repertoire in young adults before and after seasonal influenza
> vaccination. *Nat.* *Med.* **22**, 1456–1464 (2016).
>
> 2\. Kiyotani, K. et al. Characterization of the B-cell receptor
> repertoires in peanut allergic subjects undergoing oral immunotherapy.
> *J.* *Hum.* *Genet.* **63**, 239–248 (2018).
>
> 3\. Liu, S. et al. Direct measurement of B-cell receptor repertoire’s
> composition and variation in systemic lupus erythematosus. *Genes*
> *Immun.* **18**, 22–27 (2017).
>
> 4\. Kurtz, D. M. et al. Noninvasive monitoring of diffuse large B-cell
> lymphoma by immunoglobulin high-throughput sequencing. *Blood*
> **125**, 3679–3687 (2015).
>
> 5\. Riaz, N. et al. Tumor and microenvironment evolution during
> immunotherapy with nivolumab. *Cell* **171**, 934–949 (2017).
>
> 6\. Li, B. et al. Landscape of tumor-infiltrating T cell repertoire of
> human cancers. *Nat.* *Genet.* **48**, 725–732 (2016).
>
> 7\. Li, B. et al. Ultrasensitive detection of TCR hypervariable-region
> sequences in solid-tissue RNA-seq data. *Nat.* *Genet.* **49**,
> 482–483 (2017).
>
> 8\. Hu, X. et al. Landscape of B cell immunity and related immune
> evasion in human cancers. *Nat.* *Genet.* **51**, 560–567 (2019).
>
> **629**

Brief CommuniCation

> 9\. Cao, Y. et al. Potent neutralizing antibodies against SARS-CoV-2
> identified by high-throughput single-cell sequencing of convalescent
> patients’ B cells. *Cell* **182**, 73–84 (2020).
>
> 10\. Mose, L. E. et al. Assembly-based inference of B-cell receptor
> repertoires from short read RNA sequencing data with V’DJer.
> *Bioinformatics* **32**, 3729–3734 (2016).
>
> 11\. Bolotin, D. A. et al. Antigen receptor repertoire profiling from
> RNA-seq data. *Nat.* *Biotechnol.* **35**, 908–911 (2017).
>
> 12\. Chen, S.-Y., Liu, C.-J., Zhang, Q. & Guo, A.-Y. An ultrasensitive
> T-cell receptor detection method for TCR-seq and RNA-seq data.
> *Bioinformatics* **36**, 4255–4262 (2020).

13\. Mandric, I. et al. Profiling immunoglobulin repertoires across
multiple human tissues using RNA sequencing. *Nat.* *Commun.* **11**,
3126 (2020).

> 14\. Sulea, T. et al. Structure-based engineering of pH-dependent
> antibody binding for selective targeting of solid-tumor
> microenvironment. *mAbs* **12**, 1682866 (2020).
>
> 15\. Chi, X. et al. A neutralizing human antibody binds to the
> N-terminal domain of the Spike protein of SARS-CoV-2. *Science*
> **369**, 650–655 (2020).
>
> **630**
>
> **Nature** **Methods**
>
> 16\. Upadhyay, A. A. et al. BALDR: a computational pipeline for paired
> heavy and light chain immunoglobulin reconstruction in single-cell
> RNA-seq data. *Genome* *Med.* **10**, 20 (2018).
>
> 17\. Canzar, S., Neu, K. E., Tang, Q., Wilson, P. C. & Khan, A. A.
> BASIC: BCR assembly from single cells. *Bioinformatics* **33**,
> 425–427 (2017).
>
> 18\. Rizzetto, S. et al. B-cell receptor reconstruction from
> single-cell RNA-seq with VDJPuzzle. *Bioinformatics* **34**, 2846–2847
> (2018).
>
> 19\. Hagemann-Jensen, M. et al. Single-cell RNA counting at allele and
> isoform resolution using Smart-seq3. *Nat.* *Biotechnol.* **38**,
> 708–714 (2020).
>
> 20\. Zheng, G. X. Y. et al. Massively parallel digital transcriptional
> profiling of single cells. *Nat.* *Commun.* **8**, 14049 (2017).

21\. Stuart, T. et al. Comprehensive Integration of single-cell data.
*Cell* **177**, 1888–1902 (2019).

**Publisher’s** **note** Springer Nature remains neutral with regard to
jurisdictional claims in published maps and institutional affiliations.

© The Author(s), under exclusive licence to Springer Nature America,
Inc. 2021

> **NATURe** **MeTHODS**\| VOL 18 \| JUnE2021 \| 627–630 \|
> [www.nature.com/naturemethods](http://www.nature.com/naturemethods)

**Nature** **Methods**

> **Methods**
>
> **Algorithm** **overview.** TRUST4 reconstructs the immune repertoire
> in three stages: candidate reads extraction, de novo assembly and
> annotation (Fig. 1a).
>
> **Candidate** **reads** **extraction.** TRUST4 can find candidate TCR
> and BCR reads from either raw sequence files or the alignment file
> produced by aligners such as STAR22 and HISAT23. When input is an
> alignment file, if a read or its mate aligns on the V, J or C locus,
> this read is added to the candidate read set. If a read is unmapped
> and is not a candidate based on mate information, TRUST4 will test
> whether this read has a significant overlap with V, J and C genes. If
> so, this read and its mate are also candidate reads. When input is raw
> sequence files, TRUST4 applies the significant overlap criterion to
> every read or read pair to find candidate reads. To identify whether a
> read has significant overlap with one of the V, J or C genes, TRUST4
> first locates the receptor gene with the highest number of *k*-mer
> hits (default, *k*=9) from the read. TRUST4 then computes the longest
> chain from these *k*-mers to filter incompatible hits. Last, if the
> union bases of the *k*-mers in the longest chain reach threshold,
> TRUST4 will claim that the read has a significant overlap with the
> gene. The threshold is maximum(21, read_length/5+1), so data with
> shorter reads have less stringent criteria. Since TRUST4 avoids
> alignment in the candidate reads extraction algorithm, this stage is
> fast even if input data are raw sequence files.
>
> If the data have barcode information, such as 10x Genomics scRNA-seq
> data, TRUST4 also corrects the barcode, if erroneous, for each
> candidate read when given the whitelist. TRUST4 first builds barcode
> usage distribution from the first twomillion reads before correcting.
> Then, for each input barcode that is not in the whitelist, TRUST4
> finds all the neighbor barcode within one hamming distance
>
> in the whitelist (at most, fourfold barcode length) and reports the
> one that is the most frequent barcode in usage distribution. If there
> are multiple valid neighbor barcodes with the same frequency in usage
> distribution, TRUST4 will correct on the base with the lowest FASTQ
> quality.
>
> **De** **novo** **assembly.** When assembling candidate reads into
> immune receptor sequences, TRUST4 adopts the read overlap scheme.
> Cells such as plasma Bcells can generate thousands of reads for each
> recombined receptor gene, so comparison of every pair of reads to
> construct the overlap graph, as in previous versions of TRUST, is
> inefficient. TRUST4 implements a greedy extension approach by aligning
> the candidate read to existing contigs, one by one. To perform
> alignment, TRUST4 builds an index for all *k*-mers in the contigs and
> applies the seed-extension paradigm to identify alignments. TRUST4
> deems that a read overlaps with a contig if they have a highly similar
> (90% for BCR, 95% for TCR) alignment block containing at least 31-bp
> exact matches and the
>
> unaligned bases of the read are outside of the contig. Based on
> overlaps, TRUST4 will update contigs with the following rules. (1) If
> a read partially overlaps
>
> one contig, TRUST4 extends this contig; (2) if a read partially
> overlaps several contigs, TRUST4 merges corresponding contigs; and (3)
> if a read does not overlap any existing contigs, TRUST4 creates a new
> contig with this read’s sequence. When processing reads, TRUST4
> prioritizes those derived from highly expressed TCRs/BCRs. To achieve
> this, TRUST4 first counts the frequency of *k*-mers
>
> (21-mer by default) across all candidate reads. If a read comes from a
> highly expressed receptor sequence, all of its *k*-mers would be of
> high frequency in the data. Therefore, the minimum frequency of a
> read’s *k*-mers is a rough indicatation of gene abundance. TRUST4 then
> sorts the read based on the minimum *k*-mer frequency rule. The
> ordering of reads is equivalent to picking the most frequent
>
> *k*-mer as the starting point in the de Bruijn-graph-based
> transcriptome assembler, Trinity24.
>
> TRUST4 clusters reads with somatic hypermutations into the same contig
> by representing a contig as the consensus of assembled reads. Each
> position in the contig records four weights according to the number of
> reads with the corresponding nucleotide for that position. The
> consensus means that the

nucleotide for a specific position is that with the highest weight, and
the index for alignments stores the *k*-mers of the consensuses. Read
alignment takes the

> weights into account to tolerate the somatic hypermutations in BCRs.
> For example, for a particular position, if nucleotides A and T on the
> contig have high weights,
>
> it is a match if the read has nucleotide A or T. Therefore, reads with
> different somatic hypermutations can align to the same contig, which
> avoids the creation of redundant contigs.
>
> If input data are paired end, TRUST4 will use mate-pair information to
> extend the contigs. In the first round of contig assembly, due to read
> sorting and greedy extension, a contig for the abundant recombined
> gene attracts all reads from the same V, J and C genes even though
> these reads come from different recombinations. The mate-pair
> information fixes this issue by reassigning reads to the appropriate
> contigs. Reassignment will extend the contigs and update position
> weights in the affected consensus. When input data are SMART-seq,
> since there is no need to perfect assemblies for low-abundant
> sequences in a cell, TRUST4 can skip the extension to reduce running
> time.
>
> When input contains barcode information, TRUST4 will assign the read
> barcode to the contig when creating a new contig, and a read can align
> to the contigs only with that read’s barcode. As a result, two
> identical reads with different barcodes will change different sets of
> contigs. Furthermore, the read–contig
>
> Brief CommuniCation

overlap criterion is relaxed and requires 17- rather than 31-bp exact
matches in the alignment.

**Annotation.** TRUST4 aligns the assembled contigs to sequences from
the international ImmunoGeneTics (IMGT) database25 to identify V, J and
C genes. The IMGT database curates the sequences for V, D, J and
constant genes and is widely used to annotate BCR and TCR sequences,
such as in previous TRUST versions and MiXCR. Besides the sequences,
IMGT also annotates the start position of CDR3 in the V gene (104th
amino acids of the V gene in IMGT coordination). IMGT also defines the
end position of CDR3 as amino acid W/F in the amino acid motif W/FGxG in
the J gene. TRUST4 determines the CDR3

coordinate based on these IMGT conventions after identification of V and
J genes. If the contig is too short to identify the V gene, TRUST4
locates the CDR3 start position as amino acid C in the motif YYC by
testing all open reading frames.

In the final step of annotation, TRUST4 retrieves somatic hypermutated
CDR3s and estimates CDR3 abundances. If a read fully covers the CDR3 on
a contig and the CDR3 sequence from the read is different from the
consensus, TRUST4 will report the CDR3 from the read. If there is no
such read, TRUST4 directly reports the consensus CDR3 sequence. In
abundance estimation, if reads partially overlap with CDR3, each could
be compatible with several different complete CDR3 sequences. Therefore,
TRUST4 applies the expectation-maximization algorithm26, similar to that
in RSEM27, to distribute read counts iteratively to their compatible
CDR3s. For TCR, TRUST4 filters CDR3s with abundance \<5% of the most
abundant CDR3 from the same contig, to avoid sequencing errors.

For CDR3s that have only start or end positions determined in the
contigs, TRUST4 reports these as partial CDR3s and tries to extend
partial TCR CDR3s as in MiXCR. As an example of a missing start
position, the extendable partial CDR3 must overlap with the identified V
gene in the contig but cannot reach the start position. This scenario
could happen when the V gene is identified through mate-pair
information. TRUST4 then fills the missing sequences with germline
sequences of the V gene to complete the partial CDR3. In scRNA-seq,
TRUST4 also utilizes information across all cells to extend partial TCR
CDR3s. For two cells, A and B with the same V and J genes on both
chains, cell B can extend its partial CDR3 if B has a complete CDR3
identical to A’s corresponding complete CDR3, and B’s partial CDR3 is a
substring of A’s corresponding complete CDR3.

**Sequence** **data.** We tested TRUST4 on both in silico and real data.
In silico bulk RNA-seq data for evaluation of TCRs were generated using
scripts from MiXCR11, where repseqio
[(https://github.com/repseqio/repseqio](https://github.com/repseqio/repseqio))
and ART28 generated the simulated TRB and RNA-seq data. As a result,
each of the in silico RNA-seq samples contained 1,000read fragments
randomly derived from 1,000recombined TRBs. To evaluate BCR
reconstruction, we used six sets of lung cancer RNA-seq data and their
pairing BCR-seq data from our previous study8. iRepertoire processed the
BCR-seq data, and results were the gold standard for evaluation. For
SMART-seq evaluation we used three SMART-seq datasets from BALDR: AW1,
1620PV (AW2-AW3) and VH_CD19pos. For pseudobulk RNA-seq data we first
added a pseudomate for 1620PV single-end data with a sequence of one
nucleotide N. We then randomly selected 25million read pairs across all
the cells of these three samples to create

the pseudobulk RNA-seq. Finally, 56, 33 and 11% of the pseudobulk
RNA-seq data were derived from AW1, 1620PV and VH_CD19pos, respectively.
We applied the same procedure to generate psuedobulk samples with fewer
read pairs (12million, 6million, …, 2,500, 1,000). The 10x Genomics
scRNA-seq data and 10x V(D)J data were downloaded from the 10x Genomics
website.

**Performance** **evaluation.** All methods utilized were tested with
their default parameters without explicit clarification. BAM files as
input for TRUST4 were generated by STAR v.2.5.3a. In this study we used
MiXCR v.3.0.12, CATT with GitHub commit ID 0e7b462, TRUST3 v.3.0.3,
BALDR with GitHub commit ID e865b45 and BASIC v.1.5.1. All evaluations
in this work were at the nucleotide level: for example, the match of
CDR3s of TRUST4 and BCR-seq gold standard meant that their nucleotide
sequences were identical. TRUST4 can report both partial and complete
CDR3s, but we considered only complete CDR3s in evaluations.

In the TCR evaluation with in silico RNA-seq data, evaluation criteria
were based on scripts from MiXCR’s manuscript11. We added read length
150bp and ran MiXCR with default parameters. In MiXCR’s original
manuscript, the authors used the option ‘--badQualityThreshold 0' for
higher sensitivity (MiXCR_0),

and TRUST4 still found about 8% more CDR3s than MiXCR_0 on average
(Supplementary Fig. 1a). Furthermore, TRUST4 with input from FASTQ and
BAM files showed almost identical results, which demonstrated the
efficiency of the candidate extraction method. TRUST4 was also the most,
or among the most, sensitive method in assembly of CDR3s for TRB chains
with varying numbers of reads (Supplementary Fig. 1b).

For bulk RNA-seq data we mainly evaluated the performance of
reconstructing BCR heavy chains, including V, J and C gene assignments
and CDR3 sequences. We considered gene assignments in addition to CDR3
sequences in the evaluation because IGHs had different C genes as
isotypes, such as IgM, IgG1 and IgA1, and were critical in determining
antibody functions. Since CATT could not report

the C gene and TRUST3 focused only on CDR3 assembly, we omitted CATT and
TRUST3 in this evaluation. For the match of V and J genes we ignored
allele

> **NATURe** **MeTHODS**\|
> [www.nature.com/naturemethods](http://www.nature.com/naturemethods)

Brief CommuniCation

> ID. For example, if the V gene was annotated as *IGHV1-18\*01* we
> regarded it as *IGHV1-18*. In the evaluation, we excluded assemblies
> missing the V, J or C gene from TRUST4 and MiXCR. The IGH abundances
> reported by TRUST4 had a better correlation with the corresponding
> abundances in BCR-seq than MiXCR (Pearson *r*=0.57 versus 0.53 on
> average; Supplementary Fig. 2a). We further checked the
> precision–recall curve by ranking inferred IGHs by abundance (top 100,
> 500, 1,000, …), and TRUST4 consistently outperformed MiXCR across
> different thresholds (Supplementary Fig. 2a). On these real data,
> MiXCR_0 did not outperform MiXCR as in the in silico data, suggesting
> that the parameter is not effective with real data. TRUST4 with FASTQ
> and BAM input still showed identical performance across six samples in
> this real-data evaluation. We further evaluated the performance on
> CDR3 sequences only, which included results from TRUST3 and CATT.
> TRUST4 showed the highest sensitivity consistently across all six
> samples, and reported 11% more correct CDR3s than MiXCR_0, the second
> most sensitive method, with similar precision on average.
>
> The evaluations with SMART-seq data focused on whether the methods
> could reconstruct all nucleotides in the variable domain. If the
> assembled V and J sequences were shorter than gene lengths in the IMGT
> database, we regarded that as unreconstructed. The match of V or J
> sequences means that nucleotide bases were the same for regions
> annotated as V or J genes. In other words, we ignored
>
> bases before the V or after the J gene. In addition to the pseudobulk
> RNA-seq data from the three samples, we ran TRUST4 on original
> cell-level data and compared it with BALDR and BASIC on all three
> samples. We selected the top abundant heavy chain and light chain from
> TRUST4, and these were identical to either BALDR and BASIC on 272 out
> of 274chains (Supplementary Fig. 3). The comparison result indicated
> that TRUST4 can effectively reconstruct the immune repertoire from
> SMART-seq scRNA-seq data.
>
> For evaluation with 10x Genomics data, we used TCR library and IG
> library results from 10x Genomics Immune profiling (10x V(D)J) as the
> gold standard. Since the computational software CellRanger_VDJ can
> report multiple CDR3s for a cell, we regarded the most abundant CDR3
> as the true CDR3 for a chain, and the less abundant CDR3s as
> secondary. TRUST4 took the BAM file generated by CellRanger as input,
> which included the barcode information in the field “CB”. TRUST4 also
> took FASTQ files as input, and corrected the erroneous barcodes based
> on the whitelist provided in the CellRanger package. TRUST4 with FASTQ
> input reported almost identical results to that with BAM input
> (Supplementary Fig. 5c). Even though CellRanger_VDJ (v.3.1.0) was
> designed for 10x V(D)J data, we ran it on the 10x5′ scRNA-seq data
> using IMGT sequences as reference with eight cores. In our analysis of
> full-length assemblies, somatic hypermutation rate was represented by
> the proportion of matched bases (similarity) between the
>
> assembled V genes and germline sequences (Fig. 2c). When there are
> many somatic hypermutations, the similarity will be low. Besides the
> 5′scRNA-seq data, we also evaluated TRUST4 on the 3′10x Genomics PBMC
> data, with only 335cells having reconstructed CDR3s (Supplementary
> Fig. 7). We used LM22 marker genes from CIBERSORT29 to determine cell
> types.
>
> *Application* *of* *TRUST4* *on* *TCGA* *COAD* *RNA-seq* *samples.* We
> explored immune repertoire features on 466COAD RNA-seq samples in TCGA
> cohorts. To reduce the effects of somatic hypermutated CDR3s, we first
> clustered highly similar CDR3
>
> nucleotide sequences of the same length and with the same V and J gene
> assignments reported from TRUST4. We selected the similarity cutoff as
> 0.8 by comparison of similarity distribution among pairs of CDR3s
> within (intra-patient) and between (inter-patient) samples, where
> inter-patient distribution can be regarded as background random CDR3
> pair similarity (Supplementary Fig. 4a). Therefore, we defined the
> clonotype for TCR as the CDR3 sequence and that for BCR as the cluster
> with the same V and J gene assignments and similar CDR3 sequences.
> Although TRB and IGH clonalities were positively correlated with their
> respective expression (Spearman *r*=0.346 for TRB, *r*=0.085 for IGH),
> they contained additional information on TCR and BCR clonal expansion
> (Supplementary Fig. 4b). The expression for a chain is computed by the
> sum of transcripts per million (TPM) obtained from TCGA cohorts on the
> constant genes of a chain. We defined clonality as 1–(normalized
> Shannon entropy) based on the clonotype definition above.
>
> We identified that IgA1 antibody clonal expansion was related to
> patient survival in COAD. Unlike in melanoma, where IgG1 and IgA
> expression and abundance fractions were respectively positively and
> negatively associated with survival time11, we did not observe such
> association of survival time in COAD (Supplementary Fig. 4c). However,
> higher clonality of IgA1 Bcells was correlated with significantly
> shorter survival time (*P*=8.1×10–5, hazard ratio=9.14
>
> by Cox proportional hazards regression corrected by age), supporting
> the immunosuppressive property of IgA antibodies30. We hypothesize
> that the clonal expansion of IgA1 Bcells could be related to gut
> microbiota31, and future work is needed to elucidate the mechanisms
> involved.
>
> **Reporting** **Summary.** Further information on research design is
> available in the Nature Research Reporting Summary linked to this
> article.
>
> **Nature** **Methods**

**Data** **availability**

The original scripts for generation and evaluation of in silico RNA-seq
data are available at
<https://github.com/milaboratory/mixcr-rna-seq-paper>.

The six bulk RNA-seq samples for BCR evaluation are available in the SRA
repository, accession code
[PRJNA492301](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA492301),
and their matched iRepertoire data are available at
<https://bitbucket.org/liulab/ng-bcr-validate/src/master/iRep>.
SMART-seq data are available in the SRA repository, accession code
[SRP126429](https://www.ncbi.nlm.nih.gov/sra/?term=SRP126429). 10x
Genomics scRNA-seq data are available at
[https://support.10xgenomics.](https://support.10xgenomics.com/single-cell-vdj/datasets/3.1.0/vdj_nextgem_hs_pbmc3)
[com/single-cell-vdj/datasets/3.1.0/vdj_nextgem_hs_pbmc3](https://support.10xgenomics.com/single-cell-vdj/datasets/3.1.0/vdj_nextgem_hs_pbmc3),
[https://support.](https://support.10xgenomics.com/single-cell-vdj/datasets/2.2.0/vdj_v1_hs_nsclc_5gex)
[10xgenomics.com/single-cell-vdj/datasets/2.2.0/vdj_v1_hs_nsclc_5gex](https://support.10xgenomics.com/single-cell-vdj/datasets/2.2.0/vdj_v1_hs_nsclc_5gex)

and
[https://support.10xgenomics.com/single-cell-gene-expression/](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.1.0/5k_pbmc_protein_v3_nextgem)
[datasets/3.1.0/5k_pbmc_protein_v3_nextgem](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.1.0/5k_pbmc_protein_v3_nextgem).

**Code** **availability**

TRUST4 source code is available at
<https://github.com/liulab-dfci/TRUST4>. Evaluation code for this work
is available at
[https://github.com/liulab-dfci/](https://github.com/liulab-dfci/TRUST4_manuscript_evaluation)
[TRUST4_manuscript_evaluation](https://github.com/liulab-dfci/TRUST4_manuscript_evaluation).

**References**

> 22\. Dobin, A. et al. STAR: ultrafast universal RNA-seq aligner.
> *Bioinformatics* **29**, 15–21 (2013).
>
> 23\. Kim, D., Langmead, B. & Salzberg, S. L. HISAT: a fast spliced
> aligner with low memory requirements. *Nat.* *Methods* **12**, 357–360
> (2015).
>
> 24\. Grabherr, M. G. et al. Full-length transcriptome assembly from
> RNA-seq data without a reference genome. *Nat.* *Biotechnol.* **29**,
> 644–652 (2011).
>
> 25\. Lefranc, M.-P. IMGT, the international ImMunoGeneTics information
> system. *Cold* *Spring* *Harb.* *Protoc.* **2011**, 595–603 (2011).
>
> 26\. Dempster, A. P., Laird, N. M. & Rubin, D. B. Maximum likelihood
> from incomplete data via the EM algorithm. *J.* *R.* *Stat.* *Soc.*
> *Series* *B* *Stat.* *Methodol.* **39**, 1–22 (1977).
>
> 27\. Li, B. & Dewey, C. N. RSEM: accurate transcript quantification
> from
>
> RNA-seq data with or without a reference genome. *BMC*
> *Bioinformatics* **12**, 323 (2011).
>
> 28\. Huang, W., Li, L., Myers, J. R. & Marth, G. T. ART: a
> next-generation sequencing read simulator. *Bioinformatics* **28**,
> 593–594 (2012).
>
> 29\. Newman, A. M. et al. Robust enumeration of cell subsets from
> tissue expression profiles. *Nat.* *Methods* **12**, 453–457 (2015).
>
> 30\. Sharonov, G. V., Serebrovskaya, E. O., Yuzhakova, D. V.,
> Britanova, O. V. & Chudakov, D. M. B cells, plasma cells and antibody
> repertoires in the tumour microenvironment. *Nat.* *Rev.* *Immunol.*
> **20**, 294–307 (2020).
>
> 31\. Bunker, J. J. & Bendelac, A. IgA responses to microbiota.
> *Immunity* **49**, 211–224 (2018).

**Acknowledgements**

We thank B. Li and C. Wang for the helpful discussions. We also
acknowledge funding from NIH (grant U01CA226196) and China Scholarship
Council (Z.O. and Y.C.) to support this work. The study used data
generated by the TCGA Research Network that are not otherwise cited:
<https://www.cancer.gov/tcga>

**Author** **contributions**

L.S., X.H. and X.S.L conceived the project. L.S. designed and
implemented the methods. L.S., D.C., Z.O., Y.C., X.H. and X.S.L.
evaluated the methods and wrote the manuscript. All authors read and
approved the final manuscript.

**Competing** **interests**

X.S.L. is a cofounder, scientific advisory board (SAB) member and
consultant of GV20 Oncotherapy and its subsidiaries, SAB memner of
3DMedCare, consultant for

Genentech, stockholder of AMGN, JNJ, MRK and PFE and receives sponsored
research funding from Takeda and Sanofi. X.H. conducted the work while a
postdoctorate fellow at DFCI, and is currently a full-time employee of
GV20 Oncotherapy.

**Additional** **information**

**Supplementary** **information** The online version contains
supplementary material available at
<https://doi.org/10.1038/s41592-021-01142-2>.

**Correspondence** **and** **requests** **for** **materials** should be
addressed to X.S.L.

**Peer** **review** **information** *Nature* *Methods* thanks Aly Azeem
Khan, Gur Yaari and the other, anonymous reviewer(s) for their
contribution to the peer review of this work. Lin Tang was the primary
editor on this article and managed its editorial process and peer review
in collaboration with the rest of the editorial team.

**Reprints** **and** **permissions** **information** is available at
[www.nature.com/reprints](http://www.nature.com/reprints).

> **NATURe** **MeTHODS**\|
> [www.nature.com/naturemethods](http://www.nature.com/naturemethods)

<img src="./qhgfh23k.png"
style="width:3.045in;height:0.54333in" /><img src="./nfzxahel.png" style="height:0.37847in" /><img src="./f5gjeymk.png"
style="width:0.11024in;height:0.48958in" /><img src="./kks4u4yp.png"
style="width:0.13281in;height:0.54861in" /><img src="./j0wnevsd.png"
style="width:0.10764in;height:0.56163in" /><img src="./kq1pydck.png" style="height:0.42101in" /><img src="./4alnwrxj.png"
style="width:0.72656in;height:0.10764in" /><img src="./an452evd.png"
style="width:0.47743in;height:0.10851in" /><img src="./hkgyqk4v.png" style="width:0.29514in" /><img src="./hany0ega.png" style="width:0.17187in" /><img src="./f4vma5xh.png"
style="width:2.09333in;height:0.11233in" /><img src="./c0uyk5ph.png" style="width:0.21181in" /><img src="./obi0p1f4.png" style="width:0.18316in" /><img src="./qfpakosh.png"
style="width:0.40191in;height:0.10851in" /><img src="./1auz3t52.png"
style="width:0.1059in;height:0.10851in" /><img src="./nxzo2ooz.png"
style="width:0.47483in;height:0.10851in" /><img src="./rv0k4swl.png"
style="width:0.1684in;height:0.18142in" /><img src="./ymhqmdug.png"
style="width:0.11458in;height:0.19965in" /><img src="./nvxad3wq.png"
style="width:0.19184in;height:0.14497in" /><img src="./bgbznwmn.png"
style="width:0.1033in;height:0.1467in" /><img src="./aqxrc1da.png"
style="width:7.48667in;height:0.22062in" /><img src="./trbb1exq.png" style="width:0.3316in" /><img src="./u1gtonu3.png" style="width:0.43663in" /><img src="./yxwrzuzj.png" style="width:0.32812in" /><img src="./bebwwisi.png" style="width:0.10069in" /><img src="./i5kixyha.png"
style="width:0.40104in;height:0.1033in" /><img src="./34rmt2fb.png" style="width:0.15799in" /><img src="./1gcje2tw.png"
style="width:0.71962in;height:0.10851in" /><img src="./aa4qtmdo.png" style="width:0.28559in" /><img src="./i2ig3lhc.png" style="width:0.24306in" /><img src="./o2u3ysb1.png" style="width:0.19965in" /><img src="./yx1unwbz.png"
style="width:0.37413in;height:0.10851in" /><img src="./prvam4wx.png" style="width:0.19271in" /><img src="./uygwmaxa.png" style="width:0.23177in" /><img src="./ojzyzl20.png"
style="width:0.41319in;height:0.10851in" /><img src="./pqf4iqc1.png" style="width:0.45052in" /><img src="./eob3roke.png" style="width:0.13976in" /><img src="./ajhsnufb.png"
style="width:0.57552in;height:0.1033in" /><img src="./2xbhvum5.png" style="width:0.17187in" /><img src="./nqqpgu2j.png" style="width:0.65278in" /><img src="./e13toqpy.png"
style="width:0.48003in;height:0.1033in" /><img src="./bnvyewdj.png" style="width:0.15191in" /><img src="./o5w0uwfd.png" style="width:0.34983in" /><img src="./0vbjdvr2.png" style="width:0.57726in" /><img src="./mcmwjtvn.png" style="width:0.3316in" /><img src="./hcs2cref.png" style="width:0.43663in" /><img src="./13z5wdr4.png"
style="width:0.38802in;height:0.10851in" /><img src="./l2hfq1lt.png" style="width:0.82in" /><img src="./p5b1zedx.png" style="width:0.17187in" /><img src="./vuefk2vb.png" style="width:0.16059in" /><img src="./ys3x0wwa.png"
style="width:1.23229in;height:0.10851in" /><img src="./d15mwdzq.png"
style="width:0.67535in;height:0.12674in" /><img src="./24jv5wis.png" style="width:0.15191in" /><img src="./pqurvnlx.png" /><img src="./ylih5dlq.png" style="width:0.45833in" /><img src="./dpl2yfmt.png"
style="width:0.4401in;height:0.10764in" /><img src="./32c0wou3.png" style="width:0.37326in" /><img src="./v42yt4ov.png" style="width:0.19965in" /><img src="./32ghvpjh.png" style="width:0.15972in" /><img src="./g4blxyfw.png"
style="width:0.45399in;height:0.10937in" /><img src="./3ufgm5fg.png" style="width:0.26128in" /><img src="./s34u1oo1.png" style="width:0.36979in" /><img src="./3wgxupob.png" style="width:0.15885in" /><img src="./qevcxhkp.png"
style="width:0.28299in;height:0.10937in" /><img src="./2oyjlzps.png"
style="width:0.35156in;height:0.10764in" /><img src="./ry2tidtv.png" style="width:0.24653in" /><img src="./xzehw0gu.png"
style="width:0.35069in;height:0.10764in" /><img src="./azaydbab.png" style="width:0.23003in" /><img src="./tr3r50lx.png" style="width:0.21701in" /><img src="./fz2g42cv.png" style="width:0.43316in" /><img src="./htdjujnd.png" style="width:0.37674in" /><img src="./5srh10k4.png"
style="width:0.14931in;height:0.10677in" /><img src="./vlqg1jza.png" style="width:0.51302in" /><img src="./ehua52uj.png" style="width:0.17882in" /><img src="./04jl1zfa.png" style="width:0.25694in" /><img src="./lniolhxa.png"
style="width:0.34809in;height:0.10851in" /><img src="./pyqmyi3x.png" style="width:0.17535in" /><img src="./asxvwqwq.png"
style="width:0.11806in;height:0.10851in" /><img src="./hlcah0uu.png" style="width:0.13802in" /><img src="./lh34xeij.png" style="width:0.22222in" /><img src="./u0a3ijjs.png"
style="width:0.6467in;height:0.10851in" /><img src="./v0ekn0ei.png"
style="width:0.83507in;height:0.11285in" /><img src="./ohsd3hlh.png"
style="width:0.25521in;height:0.1033in" /><img src="./x3y5ixzr.png" style="width:0.38889in" /><img src="./jykqb1oj.png" style="width:0.38628in" /><img src="./2cukojcx.png" style="width:0.17187in" /><img src="./ua3f55dm.png" style="width:0.18316in" /><img src="./ny2hmhvk.png" /><img src="./s4012ipb.png" style="width:0.68924in" /><img src="./umecconj.png" style="width:0.50781in" /><img src="./zss2tm1t.png" style="width:0.4184in" /><img src="./p4a4lqqf.png" style="width:0.73611in" /><img src="./t5pftvih.png" style="width:0.26997in" /><img src="./rjc315vk.png" style="width:0.23003in" /><img src="./5hhzzddq.png" style="width:0.35851in" /><img src="./ppudykpa.png"
style="width:0.3967in;height:0.10851in" /><img src="./e1zedefn.png" style="width:0.4184in" /><img src="./wlhsxsrz.png" style="width:0.15799in" /><img src="./utmsw2ma.png"
style="width:0.34809in;height:0.10851in" /><img src="./hnejh3tf.png" style="width:0.48177in" /><img src="./nfcmtcnh.png"
style="width:0.52865in;height:0.10851in" /><img src="./ppea2h02.png" style="width:0.17882in" /><img src="./zcwb5ubv.png" style="width:0.45486in" /><img src="./r5f52g3v.png"
style="width:0.30295in;height:0.10851in" /><img src="./didnylme.png" style="width:0.22049in" /><img src="./agtj10vp.png" style="width:0.21788in" /><img src="./lpzw1dhk.png" style="width:0.4184in" /><img src="./rizws5ru.png"
style="width:0.21788in;height:0.10851in" /><img src="./mnodddke.png" style="width:0.48872in" /><img src="./3hgaqgv2.png" style="width:0.19531in" /><img src="./xc1xh4v3.png" style="width:0.52344in" /><img src="./5vnt5rlp.png" style="width:0.10156in" /><img src="./yirdhsp1.png" style="width:0.69444in" /><img src="./spc2pdg3.png" /><img src="./bb44qjk5.png" style="width:0.36111in" /><img src="./r42ozg5c.png" style="width:0.36632in" /><img src="./ubd13peh.png" style="width:0.46962in" /><img src="./rl2kqfyt.png" style="width:0.13542in" /><img src="./n3ijp4qg.png" style="width:0.74479in" /><img src="./ggbjfgko.png"
style="width:0.54687in;height:0.10851in" /><img src="./l0gcoxvk.png" /><img src="./r5i1dx2z.png" /><img src="./jyzj1w45.png" style="width:0.5in" /><img src="./m54y2nfk.png" style="width:0.30816in" /><img src="./vlywmkl0.png"
style="width:0.54687in;height:0.10851in" /><img src="./qbljiwrc.png" /><img src="./0t01kwbc.png" style="width:0.17014in" /><img src="./cfcq4ynd.png"
style="width:0.61806in;height:0.1033in" /><img src="./g2gojlva.png" style="width:0.5816in" /><img src="./jtx51poj.png" style="width:0.21615in" /><img src="./x5a5efmz.png" style="width:0.23177in" /><img src="./tv3phb0m.png" /><img src="./r3d22m4h.png"
style="width:0.46875in;height:0.10851in" /><img src="./ub3bfeh4.png" style="width:0.17187in" /><img src="./0c3jgihj.png"
style="width:0.56337in;height:0.10851in" /><img src="./0gdsiddu.png" style="width:0.13976in" /><img src="./ato1jvrk.png"
style="width:0.39844in;height:0.10851in" /><img src="./wnm1bzta.png" style="width:0.625in;height:0.1033in" /><img src="./mahhn2ex.png" style="width:0.14323in" /><img src="./m425k3ky.png"
style="width:0.54427in;height:0.10851in" /><img src="./sjvtdimz.png" /><img src="./f05h1ybu.png" style="width:0.15885in" /><img src="./rzc2cuad.png" style="width:0.45486in" /><img src="./ipkpc0t5.png" style="width:0.56337in" /><img src="./nn1mhmxf.png"
style="width:0.4401in;height:0.10851in" /><img src="./ipweg0mh.png" style="width:0.33333in" /><img src="./mftlqsoq.png"
style="width:0.46181in;height:0.10851in" /><img src="./t3vpbjyd.png"
style="width:0.19965in;height:0.11024in" /><img src="./n0n5hvfo.png"
style="width:0.34809in;height:0.10851in" /><img src="./v0p4jwag.png" style="width:0.26649in" /><img src="./z30wdgy3.png" style="width:0.23958in" /><img src="./cqmx2rok.png" style="width:0.47396in" /><img src="./105l5rux.png"
style="width:0.19965in;height:0.11024in" /><img src="./wdkwrove.png"
style="width:0.50087in;height:0.1033in" /><img src="./4fgwda1t.png"
style="width:0.55208in;height:0.10851in" /><img src="./s0jmh33f.png" style="width:0.21875in" /><img src="./bxujfkx5.png" style="width:0.43056in" /><img src="./wbef3ymx.png"
style="width:0.20052in;height:0.11024in" /><img src="./1fvo420g.png" style="width:0.42969in" /><img src="./amhzhiiz.png"
style="width:0.48611in;height:0.10851in" /><img src="./2csdspaw.png" style="width:0.51128in" /><img src="./neci4uqi.png" style="width:0.47309in" /><img src="./g3yixxz2.png" /><img src="./d3ixpxgg.png"
style="width:0.56337in;height:0.1033in" /><img src="./sk0sdaec.png"
style="width:0.19965in;height:0.11024in" /><img src="./4hqcnqte.png" style="width:0.54167in" /><img src="./rc3ppa3x.png"
style="width:0.44531in;height:0.10851in" /><img src="./p3ywjlrn.png" style="width:0.15191in" /><img src="./5kpyuzit.png" style="width:0.16319in" /><img src="./jenhwgrk.png"
style="width:0.52517in;height:0.10764in" /><img src="./5clvsy0v.png"
style="width:0.36372in;height:0.1033in" /><img src="./2d5u5jcv.png" style="width:0.16059in" /><img src="./zwkuv4hs.png" style="width:0.18663in" /><img src="./qynbzqxf.png" style="width:0.38021in" /><img src="./e0afr0qb.png"
style="width:0.19965in;height:0.10937in" /><img src="./jdtmltqi.png" style="width:0.20833in" /><img src="./i0n5o3i4.png" style="width:0.53993in" /><img src="./oixbqb3h.png"
style="width:0.44097in;height:0.1033in" /><img src="./4nb40v1c.png" style="width:0.28212in" /><img src="./0gcxcbsd.png" style="width:0.25174in" /><img src="./5ogvg4mp.png"
style="width:0.38628in;height:0.10764in" /><img src="./t04qecha.png" /><img src="./n23lxnph.png" style="width:0.42014in" /><img src="./dtr2mbis.png" style="width:0.17187in" /><img src="./024pmrsd.png" style="width:0.26042in" /><img src="./4feaomwt.png" style="width:0.27778in" /><img src="./4yz2nalv.png" style="width:0.18663in" /><img src="./kcn103ye.png" style="width:0.26997in" /><img src="./5vlt1ouf.png" style="width:0.26997in" /><img src="./myqmjxi4.png" style="width:0.81337in" /><img src="./rupyega3.png" style="width:0.15191in" /><img src="./nn2fy3ks.png"
style="width:0.41927in;height:0.1033in" /><img src="./mpdbsemz.png"
style="width:0.40538in;height:0.10851in" /><img src="./5ixu5tbc.png" style="width:0.57639in" /><img src="./455jx3ch.png" style="width:0.15885in" /><img src="./zthncmq1.png" style="width:0.31337in" /><img src="./kddqgkqq.png" /><img src="./uuayrdgn.png"
style="width:0.27691in;height:0.1033in" /><img src="./u2hedn0l.png" style="width:0.17187in" /><img src="./dg5avngq.png" style="width:0.36806in" /><img src="./ctswbeht.png" style="width:0.25347in" /><img src="./p5h4z0bh.png" style="width:0.31944in" /><img src="./5db4tfog.png" style="width:0.25in" /><img src="./zkfyevzi.png"
style="width:0.38021in;height:0.1033in" /><img src="./vusyelqq.png" style="width:0.15191in" /><img src="./0hg0ykcn.png" style="width:0.56337in" /><img src="./mdmgqrzm.png" style="width:0.17361in" /><img src="./0zouglll.png"
style="width:0.41319in;height:0.10764in" /><img src="./lelncwrk.png"
style="width:0.38889in;height:0.10764in" /><img src="./vmr1rzgx.png" style="width:0.64497in" /><img src="./gzl23bms.png" style="width:0.28559in" /><img src="./u5cgwzra.png"
style="width:0.5816in;height:0.1033in" /><img src="./zplx433l.png" style="width:0.21094in" /><img src="./dlyiwe4x.png" style="width:0.13889in" /><img src="./fo5zikud.png" style="width:0.23177in" /><img src="./gpxmafay.png" style="width:0.17361in" /><img src="./n0bj41hp.png" style="width:0.14323in" /><img src="./bhgeg2wh.png"
style="width:0.45573in;height:0.1033in" /><img src="./30sixy2d.png" /><img src="./pdjc2t2c.png" style="width:0.48524in" /><img src="./ymig5qsm.png" style="width:0.46788in" /><img src="./ukamnghw.png" /><img src="./vhvndmpe.png" style="width:0.28212in" /><img src="./2ahg0icb.png" style="width:0.22309in" /><img src="./pivxejna.png"
style="width:0.19965in;height:0.11024in" /><img src="./wfya4psd.png" style="width:0.38542in" /><img src="./syrb5vur.png" style="height:0.1033in" /><img src="./ztb4j5gw.png" style="width:0.46267in" /><img src="./wo1tr1hp.png" style="height:0.10851in" /><img src="./ec4m2rur.png"
style="width:0.47483in;height:0.10851in" /><img src="./4sy042qe.png" style="width:0.20139in" /><img src="./lic33gvo.png"
style="width:0.21701in;height:0.10851in" /><img src="./ul4yumpd.png" style="width:0.49479in" /><img src="./x1uzsawq.png" style="width:1.32986in" /><img src="./ct0tbwdv.png" style="width:0.72049in" /><img src="./13mxltqn.png"
style="width:0.51823in;height:0.13715in" /><img src="./updobkrv.png" style="width:0.13628in" /><img src="./nmrwyswp.png"
style="width:0.16753in;height:0.13542in" /><img src="./psdw2mtg.png"
style="width:0.36372in;height:0.13542in" /><img src="./v1dfaokm.png"
style="width:0.27865in;height:0.10764in" /><img src="./phoffeqs.png" style="width:0.57465in" /><img src="./sq4tcuab.png" style="width:0.28472in" /><img src="./lf455f4m.png"
style="width:1.45667in;height:0.10937in" /><img src="./emy4x1e5.png" style="width:0.42708in" /><img src="./dvyjdiek.png" style="width:0.18316in" /><img src="./ywixed52.png" style="width:0.44792in" /><img src="./3bwzwgib.png" /><img src="./s3zmy42n.png" style="width:0.38194in" /><img src="./vkc0x4sq.png" style="width:0.27951in" /><img src="./2vsoclbl.png" style="width:0.25in" /><img src="./piiersjn.png" style="width:3.22309in" /><img src="./q1xlwzum.png"
style="width:6.4168in;height:0.32705in" /><img src="./o4fm5ikl.png" style="width:0.24306in" /><img src="./y0ora3nr.png" style="width:0.18142in" /><img src="./ywsahn2p.png" style="width:0.54167in" /><img src="./pbk0dsje.png" style="width:0.20486in" /><img src="./prnnxzy2.png" style="width:0.39149in" /><img src="./fggfp1y1.png" style="width:0.1901in" /><img src="./c11glq4v.png" style="width:0.46354in" /><img src="./lsb0fd1n.png" style="width:0.25in" /><img src="./smi1fq2l.png" style="width:0.18316in" /><img src="./aaryt4oj.png" style="width:0.32292in" /><img src="./pb2cnozp.png" style="width:0.56684in" /><img src="./1cen3wvd.png" style="width:0.21528in" /><img src="./ljotxrer.png" style="width:0.47135in" /><img src="./lyjnrbra.png" style="width:0.19618in" /><img src="./5mwiqlp2.png" style="width:0.27865in" /><img src="./32ww10tc.png" style="width:0.34375in" /><img src="./1p2ktasu.png" style="width:0.22309in" /><img src="./kjpzqmhb.png" style="width:0.18142in" /><img src="./fzf35ov1.png" style="width:0.30382in" /><img src="./oru055qb.png" style="width:0.3316in" /><img src="./uhflk0mx.png" style="width:0.40538in" /><img src="./iluin5mu.png" style="width:0.28038in" /><img src="./q3hjmosc.png" style="width:0.18142in" /><img src="./jaqp34s2.png" style="width:0.30295in" /><img src="./ev5lrvx5.png" style="width:0.3316in" /><img src="./v2tpht0c.png" style="width:0.40712in" /><img src="./yiluwotf.png" style="width:0.24653in" /><img src="./dhsnahxk.png" style="width:0.28646in" /><img src="./f2nfzhgt.png" style="width:0.28125in" /><img src="./fmt310so.png" style="width:0.28646in" /><img src="./rgnljy4j.png"
style="width:6.4168in;height:0.32747in" /><img src="./00nsbpns.png" style="width:0.28733in" /><img src="./g3jtehcm.png" style="width:0.33681in" /><img src="./ssq5bylu.png" style="width:0.28646in" /><img src="./g0vzmzo2.png" style="width:0.15538in" /><img src="./4r2cy3mr.png" style="width:0.34028in" /><img src="./kdaip5oc.png" style="width:0.49653in" /><img src="./qjtsrofk.png" style="width:0.15451in" /><img src="./pxoxhdfp.png" style="width:0.25608in" /><img src="./2o24y50b.png" style="width:0.21528in" /><img src="./eu4kbyxi.png"
style="width:0.37674in;height:0.10851in" /><img src="./ytccfanw.png" style="width:0.47135in" /><img src="./3ro4nlyr.png" style="width:0.29601in" /><img src="./oz012rxo.png" style="width:0.41146in" /><img src="./gsqfgfv0.png" style="width:0.11806in" /><img src="./gxllpwbe.png" style="width:0.3724in" /><img src="./mbvwdbgm.png" style="width:0.37674in" /><img src="./1ynjlxb5.png" style="width:0.31337in" /><img src="./zbqqxy43.png" style="width:0.40712in" /><img src="./rpec5hzb.png" style="width:0.44358in" /><img src="./3lizfgjw.png" style="width:0.3967in" /><img src="./iwmw0w3g.png" style="width:0.15799in" /><img src="./otgc34yi.png" style="width:0.32031in" /><img src="./2hpgjp5k.png" style="width:2.09896in" /><img src="./go0khevz.png" style="width:0.3316in;height:0.125in" /><img src="./bhxvxnut.png"
style="width:0.27865in;height:0.10851in" /><img src="./uswrwi4n.png" style="width:0.57465in" /><img src="./3qsi0bno.png" style="width:0.28472in" /><img src="./wshonvpu.png"
style="width:0.91333in;height:0.11024in" /><img src="./mgirmf3a.png"
style="width:0.59809in;height:0.1033in" /><img src="./oxnka0cm.png" style="width:0.23872in" /><img src="./kddqor3r.png" style="width:0.34809in" /><img src="./vidyivpi.png"
style="width:1.34858in;height:0.10851in" /><img src="./q1a2jdpq.png" style="width:0.19358in" /><img src="./ggmot2ko.png" style="width:0.50868in" /><img src="./jwj45rwo.png" style="width:0.32118in" /><img src="./pexfdr3d.png"
style="width:0.36545in;height:0.10851in" /><img src="./r52tv0di.png" style="width:0.15972in" /><img src="./sebqfyxb.png"
style="width:0.45399in;height:0.11024in" /><img src="./4vb0ph3c.png"
style="width:0.60764in;height:0.10503in" /><img src="./5naa1zmb.png" style="width:0.30903in" /><img src="./ygakksoq.png"
style="width:0.52431in;height:0.10851in" /><img src="./c5bgrjhg.png" style="width:0.42708in" /><img src="./d1tiipey.png" style="width:0.2717in" /><img src="./w2s5tsjp.png" style="width:0.29427in" /><img src="./42sgpxcw.png" style="width:0.46007in" /><img src="./l11vj1cu.png" style="width:0.17969in" /><img src="./hi0gklys.png" style="width:0.18663in" /><img src="./iivxkosl.png" style="width:0.12153in" /><img src="./we43ednf.png" style="width:0.33333in" /><img src="./luwfsm4z.png" style="width:0.38108in" /><img src="./43x2l2z0.png" style="width:0.36545in" /><img src="./n1p2ihah.png" style="width:0.11458in" /><img src="./iihbgvxk.png" style="width:0.41059in" /><img src="./nqdqv1kh.png" style="width:0.17708in" /><img src="./zzze5j3s.png" style="width:0.19792in" /><img src="./g3nhhkek.png" style="width:0.4566in" /><img src="./xeivymir.png" style="width:0.18316in" /><img src="./wat2trvs.png" style="width:0.48351in" /><img src="./qjlewcl5.png" /><img src="./i1ptnw0e.png" style="width:0.15191in" /><img src="./w4vsvjsv.png" style="width:0.49045in" /><img src="./uhm31i0y.png" style="width:0.18316in" /><img src="./5d55lkh2.png" style="width:0.46528in" /><img src="./22amyrg0.png" style="width:0.31684in" /><img src="./oq5sjzgt.png" style="width:0.28038in" /><img src="./rj505sv0.png" style="width:0.12326in" /><img src="./2ijkvtzh.png" style="width:0.47309in" /><img src="./uljoixqp.png" style="width:0.15538in" /><img src="./ob4mxjxw.png" style="width:0.45399in" /><img src="./ijxmdf0t.png" style="width:0.14149in" /><img src="./3qtfsg5l.png" style="width:0.21007in" /><img src="./fovjolwm.png" style="width:0.36545in" /><img src="./hlttkf5d.png" style="width:0.18316in" /><img src="./gyjp35ql.png" style="width:0.38194in" /><img src="./ufx51wwn.png" style="width:2.42708in" /><img src="./lb4q0pki.png" style="width:0.15972in" /><img src="./fsmiohb1.png" style="width:0.10503in" /><img src="./q1e4hc54.png" style="width:0.17795in" /><img src="./lxvaogmx.png" style="width:0.36372in" /><img src="./ybeodlbb.png" style="width:0.3533in" /><img src="./qt2dvmvg.png" style="width:0.12413in" /><img src="./1tgyycyl.png" style="width:0.16319in" /><img src="./rcunsyel.png" style="width:0.45139in" /><img src="./ixhrrcpt.png" style="width:0.38021in" /><img src="./rptjp4tq.png" style="width:0.14236in" /><img src="./xacfanmt.png" style="width:0.16493in" /><img src="./hvnpwm3q.png" style="width:0.44878in" /><img src="./o5oreytc.png" style="width:0.62847in" /><img src="./rtawvvwd.png" style="width:0.15451in" /><img src="./rwnlib1z.png" style="width:0.2066in" /><img src="./actqapnz.png" style="width:0.37847in" /><img src="./0fnx0ttu.png" style="width:0.48872in" /><img src="./acboywlg.png" style="width:0.18316in" /><img src="./bzqvrgj1.png" style="width:0.38021in" /><img src="./w5bqg2h3.png" style="width:0.33333in" /><img src="./f4wgofc2.png" style="width:2.35937in" /><img src="./rr5fdvai.png" style="width:0.49826in" /><img src="./3ex2awbq.png" style="width:0.18316in" /><img src="./rqfu3mi3.png" style="width:0.38194in" /><img src="./1znw3ztz.png" style="width:0.14149in" /><img src="./sd42g1kq.png" style="width:0.16493in" /><img src="./hq0rwljc.png" style="width:0.44618in" /><img src="./sjcdz3qd.png" style="width:0.51823in" /><img src="./ptky4t51.png" style="width:0.15538in" /><img src="./ze3fq0xg.png" style="width:0.43403in" /><img src="./50sa1ys4.png" style="width:0.46007in" /><img src="./ladbmnqt.png" style="width:0.18316in" /><img src="./pudn5enh.png" style="width:0.38021in" /><img src="./ltiahrdz.png" style="width:3.99479in" /><img src="./b1ohgpuh.png" style="width:0.33333in" /><img src="./rldbcojx.png"
style="width:7.40181in;height:0.85722in" /><img src="./rnqnk0lw.png" style="width:0.15451in" /><img src="./osr3b2wo.png" style="width:2.80208in" />

<img src="./iblazp32.png" style="height:0.37847in" /><img src="./molm4un0.png"
style="width:0.11024in;height:0.48958in" /><img src="./2reqptx3.png"
style="width:0.13281in;height:0.54861in" /><img src="./hqn3y2pr.png"
style="width:0.10764in;height:0.56163in" /><img src="./g1000eub.png" style="height:0.42101in" /><img src="./cnehtbhi.png" style="width:0.40799in" /><img src="./vqlm5zez.png" style="width:0.20486in" /><img src="./j4p5c2ft.png" style="width:0.38021in" /><img src="./3blphp41.png" style="width:1.72049in" /><img src="./pqaukvii.png"
style="width:7.40181in;height:0.40514in" /><img src="./tkzv3xku.png" style="width:0.25174in" /><img src="./ef1p2upm.png" style="width:0.15538in" /><img src="./atnpvobp.png" style="width:0.21528in" /><img src="./dl5stlxc.png" style="width:0.38021in" /><img src="./gcfky4vm.png" style="width:2.80035in" /><img src="./ocderuwh.png"
style="width:0.17014in;height:0.13455in" /><img src="./rhclh1cd.png"
style="width:0.16493in;height:0.19097in" /><img src="./loxm0xqc.png"
style="width:0.19358in;height:0.13455in" /><img src="./bkv2ztf0.png"
style="width:0.15451in;height:0.16493in" /><img src="./muhng1so.png"
style="width:7.48667in;height:0.20181in" /><img src="./yijnrpgw.png" style="width:0.30295in" /><img src="./522k0omj.png" style="width:0.28472in" /><img src="./d4izd25z.png" style="width:0.15885in" /><img src="./jrmhd4r1.png" style="width:0.29167in" /><img src="./dahvx0cb.png" style="width:0.20052in" /><img src="./3jdoete4.png" style="width:0.15972in" /><img src="./mu1admss.png" style="width:0.20139in" /><img src="./1aron1hr.png" /><img src="./iybluoki.png" style="width:0.14062in" /><img src="./bynzsa4c.png" style="width:0.22483in" /><img src="./uylmq2ii.png" style="width:0.44531in" /><img src="./lyp3gx1y.png"
style="width:0.26042in;height:0.11024in" /><img src="./g3qp2ni5.png" style="width:0.15712in" /><img src="./g5ws4th4.png" style="width:0.23524in" /><img src="./mlizeqhh.png" style="width:0.20833in" /><img src="./2psqre55.png" style="width:0.15799in" /><img src="./mdtc3l0p.png"
style="width:0.5816in;height:0.1033in" /><img src="./jmt0s1om.png" style="width:0.39844in" /><img src="./lvt1khxz.png" style="width:0.31858in" /><img src="./vvzc1iwl.png"
style="width:0.34809in;height:0.10851in" /><img src="./wy53v1cn.png" style="width:0.22309in" /><img src="./dltzf024.png" style="width:0.46615in" /><img src="./o2kjhkqh.png" style="width:0.16233in" /><img src="./aahpupyb.png" style="width:0.40625in" /><img src="./enugpuow.png" style="width:0.57205in" /><img src="./a3irctpw.png" /><img src="./zakdmz2v.png" style="width:0.26562in" /><img src="./2evs25sq.png" style="width:0.40799in" /><img src="./wzkn2hn5.png"
style="width:0.50868in;height:0.10851in" /><img src="./iays5uz4.png"
style="width:0.61979in;height:0.10851in" /><img src="./y5s1khe0.png" /><img src="./qzs3fj51.png" style="width:0.71788in" /><img src="./w1dpj4vz.png" style="width:0.40625in" /><img src="./vklajolu.png" style="width:0.3316in" /><img src="./ef1emnna.png" style="width:2.30035in" /><img src="./4v5ygixu.png"
style="width:0.19792in;height:0.19531in" /><img src="./nd1xp34k.png"
style="width:0.20486in;height:0.13542in" /><img src="./5ni3hxke.png"
style="width:0.33681in;height:0.13542in" /><img src="./mld0gbd0.png"
style="width:0.17882in;height:0.16493in" /><img src="./vrnlkjhn.png"
style="width:0.21788in;height:0.13542in" /><img src="./1b1idqyk.png"
style="width:7.48667in;height:0.20222in" /><img src="./dxtngeht.png" style="width:0.34462in" /><img src="./ltmjlhxa.png" style="width:0.23872in" /><img src="./skbiqlxo.png" style="width:0.38281in" /><img src="./crccjhw0.png" style="width:0.26997in" /><img src="./5nwzlhzt.png"
style="width:0.29514in;height:0.1033in" /><img src="./3nbz5v1g.png" style="width:0.26476in" /><img src="./vdxhz3ik.png" style="width:0.16059in" /><img src="./2yljupsz.png" style="width:0.49132in" /><img src="./ked4pdn5.png"
style="width:0.4401in;height:0.1033in" /><img src="./rrgx1ria.png" style="width:0.19618in" /><img src="./xvydbohg.png" style="width:0.21701in" /><img src="./pgbhrwct.png" style="width:0.45052in" /><img src="./34l35aaa.png" style="width:0.10677in" /><img src="./lzs1vvqp.png" style="width:0.15191in" /><img src="./0pyypntc.png" style="width:0.17795in" /><img src="./0hgj24q1.png" style="width:0.39062in" /><img src="./y2du2wwk.png" style="width:0.23524in" /><img src="./bzgj40ut.png" style="width:0.49826in" /><img src="./wwaitalb.png" style="width:0.42882in" /><img src="./hcv3ign0.png" style="width:0.36372in" /><img src="./gavriig0.png" style="width:0.18316in" /><img src="./qtm42pud.png" style="width:0.3941in" /><img src="./zwjqwkcj.png" style="width:0.15278in" /><img src="./2s5qxlcc.png" style="width:0.2066in" /><img src="./dqih0rny.png" style="width:0.23437in" /><img src="./tez5nmq2.png" style="width:0.15365in" /><img src="./5feta13x.png" style="width:0.43403in" /><img src="./e041h2z2.png" style="width:0.42795in" /><img src="./0okv4hvj.png"
style="width:6.4168in;height:0.32861in" /><img src="./0ksf5voz.png" style="width:0.15538in" /><img src="./o1zzwbxz.png" style="width:0.16059in" /><img src="./ug5r20c2.png" style="width:0.3533in" /><img src="./h0eagssr.png" style="width:0.20486in" /><img src="./ooqvmaaq.png" style="width:0.5599in" /><img src="./eguqnv0p.png"
style="width:0.35764in;height:0.10851in" /><img src="./2oipmyet.png" style="width:0.17448in" /><img src="./rdn3q3v1.png" style="width:0.15451in" /><img src="./5kml0dt4.png" style="width:0.43316in" /><img src="./pnzkispr.png" style="width:0.42708in" /><img src="./bgjvmr3y.png" style="width:0.20833in" /><img src="./ut0e3d54.png" style="width:0.38889in" /><img src="./4gtecaxu.png" style="width:0.14149in" /><img src="./3xcbe5i3.png" style="width:0.18403in" /><img src="./xvmct51d.png" style="width:0.20486in" /><img src="./j0leo4pe.png" style="width:0.14236in" /><img src="./mgevgjqm.png" style="width:0.1849in" /><img src="./y1s2aket.png" style="width:0.17708in" /><img src="./qr22n3j4.png" style="width:0.23698in" /><img src="./0bg5ewyw.png" style="width:0.14062in" /><img src="./rnqy311i.png" style="width:0.29514in" /><img src="./h0w5njlj.png" style="width:0.24826in" /><img src="./vhiaye0e.png" style="width:0.4184in" /><img src="./ckpdqkga.png" style="width:0.24306in" /><img src="./q52uylw0.png" style="width:0.36545in" /><img src="./xcsf0wra.png" style="width:0.3533in" /><img src="./tj0x41by.png"
style="width:6.4168in;height:0.32903in" /><img src="./kyld4ccc.png" style="width:0.16493in" /><img src="./ru2iaf4u.png" style="width:0.27517in" /><img src="./l2pjerkl.png" style="width:0.27951in" /><img src="./5tbjhbq0.png" style="width:0.21528in" /><img src="./rdo23see.png" style="width:0.50174in" /><img src="./xweh0omf.png"
style="width:6.4168in;height:0.20858in" /><img src="./r02oaeu3.png" style="width:0.35503in" /><img src="./xwth5ohm.png" style="width:0.3533in" /><img src="./bufcvkft.png" style="width:0.18142in" /><img src="./a02sgk1f.png" style="width:0.44184in" /><img src="./h0rcr2in.png" style="width:0.13976in" /><img src="./1t1wlbfv.png" style="width:0.61198in" /><img src="./ds4eh5ao.png" /><img src="./ll53gjvf.png" style="width:0.63628in" /><img src="./peavswj0.png" style="width:0.40538in" /><img src="./43iilaom.png"
style="width:0.53819in;height:0.10764in" /><img src="./2vccew0x.png"
style="width:6.4168in;height:0.20878in" /><img src="./1h2k5lmg.png" style="width:0.35503in" /><img src="./stpefpwy.png" style="width:0.3533in" /><img src="./hagn5xvy.png" style="width:0.18142in" /><img src="./l45bfhup.png" style="width:0.44184in" /><img src="./lnc4ygdz.png" style="width:0.13976in" /><img src="./jj52f4ag.png" style="width:0.61198in" /><img src="./1nvj2l3a.png" /><img src="./i3njtrip.png" style="width:0.63628in" /><img src="./j55fczgl.png" style="width:0.40538in" /><img src="./vccybpkp.png" style="width:0.73351in" /><img src="./ibimf3wa.png"
style="width:6.4168in;height:0.20899in" /><img src="./manzpr0w.png" style="width:0.35503in" /><img src="./dsksfw4p.png" style="width:0.3533in" /><img src="./144seh4w.png" /><img src="./hdt5fwev.png" style="width:0.14149in" /><img src="./hrxmh4q3.png" style="width:0.18316in" /><img src="./vcljgz0c.png" style="width:0.20399in" /><img src="./kp1no22r.png" style="width:0.25608in" /><img src="./st03r0ua.png" style="width:0.35069in" /><img src="./ccnj4imd.png"
style="width:0.38715in;height:0.10851in" /><img src="./ggmdor4y.png"
style="width:0.24392in;height:0.17969in" /><img src="./rwica5vn.png"
style="width:0.15451in;height:0.16493in" /><img src="./dmtytxfy.png"
style="width:0.20747in;height:0.19444in" /><img src="./00vh21ig.png" style="height:0.13455in" /><img src="./a1v3u4bv.png"
style="width:0.16319in;height:0.19097in" /><img src="./w3yxhy04.png"
style="width:0.17535in;height:0.13108in" /><img src="./hcyawp15.png"
style="width:0.3316in;height:0.16493in" /><img src="./czgxit5a.png" style="height:0.13455in" /><img src="./a0jq2ur2.png"
style="width:0.14149in;height:0.17101in" /><img src="./cx4f4czn.png"
style="width:7.48667in;height:0.23049in" /><img src="./cmd0kxkv.png"
style="width:0.17535in;height:0.13108in" /><img src="./kuufmjoo.png" style="height:0.13455in" /><img src="./odlsersv.png" style="height:0.13455in" /><img src="./4kxjpl11.png"
style="width:0.17535in;height:0.13108in" /><img src="./tnsx5t0p.png"
style="width:0.20486in;height:0.16493in" /><img src="./wbjf35ew.png"
style="width:0.25347in;height:0.19358in" /><img src="./z0jpkkze.png" style="height:0.13455in" /><img src="./ggmea42j.png" style="width:0.14236in" /><img src="./3jo5ovkt.png" style="width:0.31163in" /><img src="./44robffv.png" style="width:0.51302in" /><img src="./cx5wzu4l.png" style="width:0.20573in" /><img src="./ebbnwefl.png" style="width:0.33333in" /><img src="./d5dv012t.png" style="width:0.25174in" /><img src="./io1obu4f.png" style="width:0.23351in" /><img src="./qmzmv2jg.png" /><img src="./s3i1g4ki.png" style="width:0.42274in" /><img src="./0fqpyffy.png" style="width:0.57639in" /><img src="./ye310uco.png" style="width:0.34635in" /><img src="./s2zmda1s.png" style="width:0.15451in" /><img src="./bwreub4c.png" style="width:0.38194in" /><img src="./cskgprwf.png" style="width:0.19705in" /><img src="./xynlbkhs.png" style="width:0.23698in" /><img src="./0pyd5yut.png" style="width:0.32986in" /><img src="./0txptwvz.png" style="width:0.22309in" /><img src="./wns1wvpj.png" style="width:0.34115in" /><img src="./n4omb0t4.png" style="width:0.36979in" /><img src="./ihgsbjqq.png" style="width:0.19878in" /><img src="./ok1hrrnx.png" style="width:0.37847in" /><img src="./1uq25ucz.png" style="width:0.30035in" /><img src="./ak3hubuc.png" style="width:0.33681in" /><img src="./enc1mxvp.png" style="width:0.22483in" /><img src="./x4chdtzt.png" style="width:0.35503in" /><img src="./hzgnmdcu.png" style="width:0.19878in" /><img src="./a0v054qs.png" style="width:0.25868in" /><img src="./amtgtgzl.png" style="width:0.23177in" /><img src="./czztajax.png" style="width:0.11458in" /><img src="./qll5uv03.png" style="width:0.1875in" /><img src="./mrgbavgr.png" style="width:0.30469in" /><img src="./j5jdfh20.png" style="width:0.19792in" /><img src="./0unjiwmk.png" style="width:0.39497in" /><img src="./czce4tu5.png" style="width:0.18663in" /><img src="./e4jo23xq.png" style="width:0.14149in" /><img src="./hqcf2lvn.png" style="width:0.51649in" /><img src="./w2il5jvl.png" style="width:0.30642in" /><img src="./xou1j432.png" style="width:0.28125in" /><img src="./oswylie3.png" style="width:0.38628in" /><img src="./lcqplish.png" style="width:0.41667in" /><img src="./cvm5zmzm.png"
style="width:0.55816in;height:0.10677in" /><img src="./3bsnihbi.png" style="height:0.1033in" /><img src="./v1p2smcx.png"
style="width:0.79167in;height:0.13368in" /><img src="./4ef4kud5.png"
style="width:0.4783in;height:0.11632in" /><img src="./c04yagyo.png"
style="width:0.15017in;height:0.10677in" /><img src="./uzyhu3je.png" style="width:0.40365in" /><img src="./r2khqzbv.png" style="width:0.15885in" /><img src="./fnnqmeon.png"
style="width:0.26649in;height:0.10851in" /><img src="./nlqtddo0.png" style="width:0.46962in" /><img src="./5brkjju0.png" style="width:0.45399in" /><img src="./zcusguma.png" style="width:0.13325in" /><img src="./h1r1atfi.png" style="width:0.19184in" /><img src="./hue1i3dr.png" style="width:0.6224in" /><img src="./0ilglk0d.png" style="width:0.15538in" /><img src="./0j0ynbka.png" style="width:0.53646in" /><img src="./wiohyrwp.png" style="width:0.34635in" /><img src="./kickgb5i.png" style="width:0.15451in" /><img src="./nqhtetvl.png" style="width:0.23698in" /><img src="./sjpvds35.png" style="width:0.44705in" /><img src="./lmmsbm3o.png" style="width:0.30556in" /><img src="./ovhhfvzz.png" style="width:0.36806in" /><img src="./0lvkzqeb.png" style="width:0.51476in" /><img src="./jaatc2s4.png" style="width:0.2947in" /><img src="./ya1hw0lo.png" style="width:0.18316in" /><img src="./yqtawvra.png" style="width:0.18229in" /><img src="./23ic3huh.png" style="width:0.36806in" /><img src="./02k4q1se.png" /><img src="./lom1moq2.png"
style="width:0.53299in;height:0.10677in" /><img src="./vfqg1lao.png"
style="width:0.15017in;height:0.10677in" /><img src="./w5g3mqex.png" style="width:0.40278in" /><img src="./yutxshau.png" style="width:0.15799in" /><img src="./22ozr32k.png"
style="width:0.26649in;height:0.10851in" /><img src="./kdgcmmev.png" style="width:0.37674in" /><img src="./yamvlt1x.png" style="width:0.19792in" /><img src="./usazasgy.png" style="width:0.44965in" /><img src="./x4jhpfqg.png" style="width:0.46181in" /><img src="./lqzhhb20.png" style="width:0.60503in" />
