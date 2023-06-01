# Single cell multimodal omics

Reviewed: No
Type: Slides

All bulk omics are reaching a single-cell resolution, currently with variable degrees of success. Some examples of trans-omics encompass genome, transcriptome, proteome and metabolome. For instance, single-cell applies both to *epigenetic* (DNA and histone modifications) and *epitranscriptomics* (the ensemble of all RNA modifications).

Single cell multimodal omics was classified as method of the year 2019 by Nature.

## Single-cell multiomic profiling

Linking measurements from different omics layers has the potential to reveal regulatory and functional mechanisms underlying cell behaviour in healthy development and disease. In Figure …, we can observe the possible connections among omics techniques and their interplay in determining cell phenotype and functional potential as well as molecular cell identity.

![Screenshot.png](Single%20cell%20multimodal%20omics%206a27bf9c285c4cdd9cc4eea550d22dc6/Screenshot.png)

# Transcriptome & proteome

Proteins determine much of cell behaviour and mRNA expression levels are a weak proxy of protein expression (post-transcriptional and translation regulation). In addition, there is no method for protein sequence amplification (measurements are dependent on antibody-based detection or mass spectrometry).

## **DNA barcoded antibodies used in single-cell multiomics**

The use of **DNA-barcoded antibodies** is applied to convert the detection of proteins into a quantitative and sequenceable readout. Antibody-derived tags (**ADT**) act as **synthetic transcripts,** captured during oligodT-based scRNA- seq library preparation protocols (e.g. 10x, Drop-seq). Cells labelled with panels of antibodies, each tagged with a specific ADT, **captured in parallel with the mRNA from the same cell** following lysis.

## CITE-seq (Cellular Indexing of Transcriptomes and Epitopes) (2017)

CITE-seq uses DNA-barcoded antibodies to convert detection of proteins into a quantitative, sequenceable readout. Antibody-bound oligos act as synthetic transcripts that are captured during most large-scale oligo dT-based scRNA-seq library preparation protocols. 

![cite-seq.com](Single%20cell%20multimodal%20omics%206a27bf9c285c4cdd9cc4eea550d22dc6/Screen_Shot_2023-02-22_at_20-37-02.png)

cite-seq.com

**Pros:** immuno-phenotyping of cells allows a potentially limitless number of markers unbiased transcriptome analysis.

**Cons:** specific antibodies necessary, for surface protein markers

Transcriptome-based clustering of 8,005 CITE-seq single-cell expression profiles of Cord Blood Mononuclear Cells (library of 13 markers) ADT: antibody-derived tags.

## 2018: Cell Hashing with barcoded antibodies

Oligo-tagged antibodies against **ubiquitously expressed** surface proteins.Label cells from distinct samples, which can be **pooled** (cost reduction). **Detection of doublets** between samples (or within a sample, if used without pooling).

Other techniques based on barcoded antibodies:

- 2017: **REAP-seq** (RNA Expression and Protein Sequencing assay). Peterson et al. Nature biotechnology (library of 82 barcoded antibodies)
- 2017: **Abseq** (Ultrahigh- throughput single cell protein profiling with droplet microfluidic barcoding). Shahi et al, Scientific reports (example with 2 antibodies, B cells vs T cells)

## ****Single-cell mass-spectrometry approaches (proteome only)****

Single cells are isolated in individual wells ****and lysed. Proteins are digested to peptides. Peptides from each single cell are **covalently labeled** (barcoded) with isobaric tandem-mass-tags (TMT). 

![*Specht et al, Genome Biology 2021*
](Single%20cell%20multimodal%20omics%206a27bf9c285c4cdd9cc4eea550d22dc6/Screenshot%201.png)

*Specht et al, Genome Biology 2021*

SCoPE2 was applied to monocytes differentiated into macrophages. Quantified over 3042 proteins in 1490 single monocytes and macrophages in ten days of instrument time.

# Transcriptome & genome

Somatic variation is the genomic diversification between cells, deviating from the “original” genome of the zygote. Somatic variation can be:

- *programmed*: during the maturation of lymphocytes, programmed rearrangements in B and T cell receptors are carried out e.g. rearrangements of the V(D)J regions in B and T lymphocytes to produce diverse and specific antibodies and T cell receptors.
- *spontaneous*: accumulations of SNVs or CNVs of chromosome regions during development and aging

Potentially pathogenic when variants confer a competitive advantage to the cell, leading to clonal expansion and the formation of malignant or cancerous clones.

## DR-seq (2015)

DR-seq (gDNA-mRNA sequencing) is based on sequencing of genomic DNA and messenger RNA from the same cell. There are two different strategies for amplifying:

- PCR for gDNA amplification
- mRNA specific second-strand synthesis (polyA selection, only the 3’ end is sequenced)

Coverage depth in a single cell should be 2, if we see more is produced by amplification - differently from RNA-seq. There is a positive correlation between copy number levels and average expression of all genes located in the region of interest e.g. chromosome 8.

![Screenshot.png](Single%20cell%20multimodal%20omics%206a27bf9c285c4cdd9cc4eea550d22dc6/Screenshot%202.png)

## G&T-seq

Parallel sequencing of single-cell genomes and transcriptomes occurs through physical separation of genomic DNA and poly(A)+ mRNA after cell lysis, followed by separate amplification, library preparation and sequencing.
Full-length mRNA sequencing can be obtained following Smart-seq2 protocol.

**Workflow:**

1. isolate cells
2. lysis
3. using oligo-dT probes isolate mRNA content and divide it in different plates with respect to DNA
    1. RNA: whole transcriptome amplifications and full length sequencing
    2. DNA: whole genome amplification and detection

The method was applied to differentiated neurons derived form trisomy 21 induced pluripotent stem cells and control disomy 21 iPSCs. By visual inspection, we can clearly see a gain in chromosome 21 with respect to control, even though some error is identified (…). The positive result is that chromosome 21 genes are more expressed with respect to control. Transcriptomics variation was also observed in other autosomes, which could be due to both biological and technical variation.

![*Macaulay et al, Nature Methods 2015*
](Single%20cell%20multimodal%20omics%206a27bf9c285c4cdd9cc4eea550d22dc6/Screenshot%203.png)

*Macaulay et al, Nature Methods 2015*

In theory these approaches can be used to detect fusions, while SNVs are harder to be found at low coverage. There could be location biases, but this is common for different amplification techniques.

## Physical separation of the nucleus and cytoplasm

In this case there is no separation of DNA and RNA, but rather of the nucleus (DNA-seq) and cytoplasm (RNA-seq). Cell separation is performed by sorting, but the core of the technique a gentle lysis and separation in two plates. Full length mRNA sequencing is comparable with SMART-seq2, while for DNA tagmentation is performed.

![*Zachariadis et al, Molecular Cell 2020*](Single%20cell%20multimodal%20omics%206a27bf9c285c4cdd9cc4eea550d22dc6/Screenshot%204.png)

*Zachariadis et al, Molecular Cell 2020*

## Genotyping of Transcriptomes (GoT) (2019)

![*Nam et al, Nature 2019*](Single%20cell%20multimodal%20omics%206a27bf9c285c4cdd9cc4eea550d22dc6/Screen_Shot_2023-02-22_at_20-38-39.png)

*Nam et al, Nature 2019*

GoT integrates *genotyping* (detection of somatic mutations gene of interest) with single cell RNA sequencing. The main advantage of this technique is the identification of malignant cells associated with a specific mutation(challenging in the absence of specific surface markers). Additionally, this method works with droplet-based methods, greatly enhancing the throughput.

We can connect the presence of a mutation with markers and cell abundance. The real advantage is the amplification of the locus of interest reducing the dropout rate.

## Adaptive immune response

Both lymphocytes express B- and T-cell receptors. The soluble part in both receptor is an antibody able to recognize antigens. The transmembrane domain is common, while the top part is a variable region. Looking at gDNA of each cell, we have a V(D)J region where random stochastic rearrangement occur, producing variability in the amino acid region of the receptor. Once one of these cells meets a pathogen, we have the activation of the response and increase in antibody targeting a specific antigen. Through the sequence of the VDJ region, we can identify the abundance of the pathogen e.g. COVID. 

![ *https://www.10xgenomics.com*
](Single%20cell%20multimodal%20omics%206a27bf9c285c4cdd9cc4eea550d22dc6/Screenshot%205.png)

 *https://www.10xgenomics.com*

### Immune profiling protocol

The most used platform (10x) has an **immune profiling protocol** implementing this technique through transcriptome 5’ end sequencing for paired T- and B-cell receptors. After barcoding, libraries are split for the targeted analysis of B/T cell receptor transcripts, (based on primers specific for the sequence of the constant region) and the analysis of all other transcripts (standard 10x 5’ sequencing).

Sequencing of TCR alpha and beta chain allows cell grouping according to shared sequences. It is possible to predict the pathogen to which the receptor will react e.g. the clonal expansion is recognizing a sequence of a particular peptide of cytomegalovirus.

# Transcriptome & epigenome

The compactness of DNA is regulated by protein modifications that act on chromatin accessibility and therefore have an effect on gene transcription. Epigenetic modifications contribute to cellular heterogeneity and regulation of gene expression during development, lineage determination and response to dynamic stimuli. In single cell epigenomics we can study:

- DNA methylation
- chromatin accessibilty

## DNA  methylation

DNA methylation was the first discovered epigenetic mark, which consists on an addition of a methyl group on cytosine residues of the dinucleotide **CpG.** Methylation is implicated in repression of transcriptional activity. **Bisulphite sequencing** is a technique based on treating DNA with sodium bisulphite, which allows the conversion of cytosines to uraciles through alkylation in the case of unmethylated cytosine - while methylated cytosine will remain unchanged. By analyzing mismatches with a reference genome sequence, we can distinguish methylated from unmethylated cytosines.

### scM&T-seq (2016)

Single-cell methylome & transcirptome is a derivation of G&T-seq, based on the physical separation of DNA and RNA:

- Smart-seq2 on mRNA
- Bisulphite sequencing on DNA

The technique was applied to mouse embryonic stem cells in the first publication of the technique. Several single cell genome and transcriptome methods can be extended to study the methylome by applying bisulphite sequencing.

![*Angermueller et al, Nature Methods 2016*
](Single%20cell%20multimodal%20omics%206a27bf9c285c4cdd9cc4eea550d22dc6/Screenshot%206.png)

*Angermueller et al, Nature Methods 2016*

## Chromatin accessibility

### ATAC-seq (2013)

ATAC-seq (Assay for Transposable Accessible Chromatin) is designed to identify open chromatin regions in the genome by high-throughput sequencing. The technique is based on the activity of a key enzyme, **Tn5 Transposase**, which performs a “cut and paste” procedure and can be engineered into a hyperactive mutant form. Tn5 transposases preferentially insert into open chromatin sites, cut and add two sequencing primers (t*agmentation*). The obtained fragments will be nucleosome free or nucleosome containing. The technique is very sensitive i.e.  works with small amounts of material and faster than alternative methods. 

![*Buenrostro et al, Nature Methods 2013*
](Single%20cell%20multimodal%20omics%206a27bf9c285c4cdd9cc4eea550d22dc6/Screenshot%207.png)

*Buenrostro et al, Nature Methods 2013*

**Key features to identify a successful ATAC-seq experiment**:

- fragment size peaks should follow a clear periodicity of ~200 bp length, corresponding to multiples of nucleosome-wrapped sequence size
- at the genome-wide level, fragments will be enriched in open chromatin regions:
    - promoter flanking regions
    - transcribed genes

### Single cell ATAC-seq (2015)

ATAC-seq was adapted to single cell resolution with a fluidigm device (not yet droplet-based, lower throughput). We can analyze a snapshot of a genome portion comparing fragment peaks generated from bulk ATAC-seq to aggregate single cell ATAC-seq peaks (obtained from 254 cells) and verify that there is a high correlation among peaks, which are located close to the TSS. It is possible to map fragments on a single cell profile, which will have a range from 0 to 2 fragments per position. The protocol can be adapted to any single cell approach e.g. droplet-based or combinatorial indexing. In 2018, 10x developed 10x single-cell ATAC-seq with a similar workflow; the main difference is the usage of nuclei instead of cells and the capture procedure  without UMI (directly designed for transposase adapters). 

“Why always chromosome 19?”

![*Buenrostro et al, Nature 2015*](Single%20cell%20multimodal%20omics%206a27bf9c285c4cdd9cc4eea550d22dc6/Screenshot%208.png)

*Buenrostro et al, Nature 2015*

### Peak calling

In order to reduce the dimensions, one of the first steps in the analysis is peak calling, which can be performed by Cell Ranger on 10x data. The number of transpositions events is evaluated in pseudo-bulk, followed by a smoothing procedure and signal threshold evaluation to discriminate actual signal from noise e.g. zero negative binomial model. The outputs of this step are:

1. indexed fragment file collecting all the reads, even if they do not map to a peak, which can be used for quality control.
2. large sparse matrix, where each row is a genomic range corresponding to a peak and the column the number of fragments per cell

The main challenges in comparison to scRNA are:

- more sparse data
- near-binary data (open vs closed, we are not really interested on how much)
- non-fixed feature set (variable features according to number of peaks per cell)
- order of magnitude more features

### Enrichment

After peak calling, it is possible to perform DNA sequence enrichment on motif to identify transcription factor binding sites. Regulatory relationship among enhancers and promoters can  be studied through co-accessibility networks and gene variants can be studied through genetic variant enrichment.

While working with Seurat, the best approach is to apply **Signac**, a Seurat extension for the analysis of ATAC-seq data. 

## SNARE-seq (2019)

Single-Nucleus chromatin Accessibility and mRNA Expression is a high-throughput sequencing of the transcriptome and chromatin accessibility in the same nucleus. In this setting, there is no need for probabilistic mapping of single-cell clusters from separate analyses. The technique was applied to droplet-based method, then adopted by 10x (Multiome profiling).

Capturing beads contain two adapters, one for ployA selection for mRNA and another complementary to Tn5. The advantage is that probes have the same barcode, allowing easy mapping of both libraries and digital counting matrix with common columns.

![Chen et al, Nature Biotechnology 2019
******](Single%20cell%20multimodal%20omics%206a27bf9c285c4cdd9cc4eea550d22dc6/Screenshot%209.png)

Chen et al, Nature Biotechnology 2019
******

In the original publication the method was applied to ~5000 nuclei from 5 mice collected from neonatal cerebral cortex. The clusters have been annotated to known groups of neuronal cells based on gene expression (reference) and then mapped to the UMAP obtained from chromatin data analysis. The obtained result is more or less consistent with the expectation and allows the study of biological implications of chromatin regulation on gene expression. Alternatively, it is possible to use both approaches simultaneously to identify clusters.

![***Chen et al, Nature Biotechnology 2019***
](Single%20cell%20multimodal%20omics%206a27bf9c285c4cdd9cc4eea550d22dc6/Screenshot%2010.png)

***Chen et al, Nature Biotechnology 2019***

10x method for simultaneous detection of gene expression and chromatin state from the same cell can aid in identifying activation of gene expression even though expression data does not show relevant patterns. In addition, putative regulatory elements directly linked to a gene of interest can be inferred thanks to the correlation of gene expression with the presence of accessible peaks in clusters.

![*https://www.10xgenomics.com*
](Single%20cell%20multimodal%20omics%206a27bf9c285c4cdd9cc4eea550d22dc6/Screenshot%2011.png)

*https://www.10xgenomics.com*

# ****Transcriptome & “CRISPR” perturbations****

CRISPR technology can be applied to inactivate regions of interest and study perturbation effects (also called *reverse genetics* approach).

## Perturb-seq (2016)

Perturb-seq combines single-cell RNA sequencing with CRISPR- based perturbations (inactivation of genes). The main aim is to map the transcriptional effects of genetic perturbations and identify regulatory circuits. By applying single cell, we can detect the expression level of each cell and identify which kind of CRISPR guide was present in each cell i.e. which gene was inactivated. 

![Dixit et al, Cell 2016](Single%20cell%20multimodal%20omics%206a27bf9c285c4cdd9cc4eea550d22dc6/Screen_Shot_2023-02-22_at_20-40-56.png)

Dixit et al, Cell 2016

The mRNAs are captured with a cell barcode (**CBC**) and matched to sgRNAs and paired guide barcode (**GBC**).

![Screenshot.png](Single%20cell%20multimodal%20omics%206a27bf9c285c4cdd9cc4eea550d22dc6/Screenshot%2012.png)

## CRISP-seq (2016)

CRISP-seq combines single-cell RNAseq with CRISPR-based perturbations (editing/inactivation of genes). It allows to discover interactions and redundancies between developmental and signaling- dependent factors. The technique exploits mathematical modeling and immune stimulation factors.

![*Jaitin et al, Cell 2016*
](Single%20cell%20multimodal%20omics%206a27bf9c285c4cdd9cc4eea550d22dc6/Screenshot%2013.png)

*Jaitin et al, Cell 2016*

### 10x CRISPR screen

In addition, 10x now provides single-cell CRISPR screens, which provide information on:

1. guide RNAs
2. whole transcriptome profiling
3. cell surface protein expresssion
4. isolation of paired immune receptor sequences

# Challenges and opportunities in single-cell multimodal omics

![*Zhu C et al. Nat Methods. 2020*](Single%20cell%20multimodal%20omics%206a27bf9c285c4cdd9cc4eea550d22dc6/Screenshot%2014.png)

*Zhu C et al. Nat Methods. 2020*

Chromatin accessibility layer is more sparse than the genomic one, so the combination of multiple layers can aid in achieving a better separation. This instead in not an issue with tagged antibodies, but on the other hand we are working with restricted information and mainly on surface proteins. Epigenetic information on histone marks and TF binding, even though not in the same cell, can be combined to study regulatory networks.