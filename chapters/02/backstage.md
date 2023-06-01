# SC3

Date: March 28, 2023
Reviewed: No
Type: Lecture

# Challenges in single cell workflow

![Screenshot.png](SC3%20b0c718fdb656445abd059bcb944703a4/Screenshot.png)

## Sources of experimental variability

1. human factor
2. cell suspension preparation
3. cDNA synthesis
4. PCR
5. NGS

## Single-cell suspension

RNA degrades very fast, so it is required to quickly transfer separated samples to the sequencing machine. A good quality single-cell suspension should have high viability, no dead or dying cells, cell debris removed and unaltered transcriptional profile.

Fresh samples can be obtained from:

- suspension cell lines e.g. blood
- adherent cell lines →  trypsin digestion is usually applied for dissociation, but it is possible to  choose among other enzymes e.g. Accutase or TrypLe. It is then necessary to resuspend cells to avoid aggregates. Wide-pore tips are used to avoid pression and stress induction, leading to an overall higher cell survival.
- dissociation of tissues

Cells are different in size, gene expression profile and are characterized by peculiar extracellular matrix composition, which affects dissociation techniques. 

### Single cell suspension from solid tissue

1. tissue dissection
2. rinse with PBS to get rid of unwanted components
3. mechanical dissociation to increase the contact area with cell and enzyme
4. add digestive enzymes (go to Tissue Atlas to choose the best option for the tissue of choice)
5. incubate considering temperature, time and use of shaker
6. remove undigested portions of tissue through a filter
7. evaluate quality of single cell preparation

This protocol is time consuming and only encompasses a small piece of the pipeline; it is possible to rely on commercial solutions, as Miltenyi Biotec automated tissue dissociator with additional dead cell removal kit (through magnetic beads) or debris removal solution. 

![STAR Protocols 2, 100989, December 17, 2021
****](SC3%20b0c718fdb656445abd059bcb944703a4/Screenshot%201.png)

STAR Protocols 2, 100989, December 17, 2021
****

## Single cell sorting

While working with rare cells or genetically modified cells, it is required to apply **FACS** (Fluorescent Activated Cell Sorting), where cells are fluorescently labelled and driven by pressure in a tube where laser excites cells and a detector reads chemical properties. Finally, cells are sorted through deflector plates, included in a droplet and inserted in collection tubes. The overall procedure is long and stress-inducing, hence cells could have changes in gene expression profile (possible solution: freeze a part of the samples for final comparison).

![G Pfister et al., J Immunol Meth 487 (2020) 112902
****](SC3%20b0c718fdb656445abd059bcb944703a4/Screenshot%202.png)

G Pfister et al., J Immunol Meth 487 (2020) 112902
****

**Microfluidic sorters** provide a commercial solution to the issue. They work in a sterile environment, are based on disposables, avoid cross contaminations and apply a gentle analysis. Through the *actuator,* cells are mechanically accompanied to the right channel without stress induction. Examples:

- LeviCell, LevitasBio: label free system exploiting the magnetic properties of nanoparticles in cells for sorting. Drawback: small input channels, only a small number of cells is allowed.

## Sample collection and preservation

With the advent of molecular pathology, hospitals started to collect fresh and frozen tissue in addition to FFPE-preserved tissue traditionally used in immunohistochemistry.  From these sample, it is possible to isolate both nuclei and cells according to the techniques. In order to maintain samples in time, freezing is a feasible solution, as cells can be cultured again after many years.

![Screenshot.png](SC3%20b0c718fdb656445abd059bcb944703a4/Screenshot%203.png)

Keep in mind that different assays require different input materials:

- protein expression: whole cells
- gene expression: whole cell nuclei
- chromatin accessibility: nuclei

Determining your required input can help you determine what collection and preservation method is most suitable for your sample type.

# 10xGenomics

Three main pillar technologies:

1. Chromium Single Cell: NGS readout
2. Visium Spatial: NGS readout
3. Xenium In Situ: microscopy-based readout

In **feature barcode technology** we aim at combining different cells in the same sequencing run through cell multiplexing, cell surface protein or CRISPR screening. Older machines run in manual mode, while the more recent Chromium Connect is automated and able to perform gene expression immune profiling.

## Feature barcode technology

In the standard approach, we find 3.6M barcodes, while in low throughput (a bit cheaper) 9.2K barcodes. The capture sequence is directed to antibody or guide RNA for gene activity disruption.

![Screenshot.png](SC3%20b0c718fdb656445abd059bcb944703a4/Screenshot%204.png)

Samples are inserted in a chip ..

![Screenshot.png](SC3%20b0c718fdb656445abd059bcb944703a4/Screenshot%205.png)

**Bio-Rad vs 10x Genomics patent war:** develoment of digital PCR by Quantalife, bought by Bio-Rad. 

## Single cell gene expression workflow

1. GEM
    
    ![Screenshot.png](SC3%20b0c718fdb656445abd059bcb944703a4/Screenshot%206.png)
    
2. GEM breakage
    
    ![Screenshot.png](SC3%20b0c718fdb656445abd059bcb944703a4/Screenshot%207.png)
    
3. cDNA amplification and library construction
    
    ![Screenshot.png](SC3%20b0c718fdb656445abd059bcb944703a4/Screenshot%208.png)
    

## Chip-based single cell separator

Chip-based single cell separator → iCELL8 cv (TAKARA) 72x72 wells, check the quality by inspection through a microscope. Very low volumes, we should have no air in the solution. Upper limit for nozzle is higher than 10x, so it is possible to analyze bigger cells e.g. Schwann cells or cardiomyocytes. With this technique, reagents are only dispensed in selected wells. 

![Screenshot.png](SC3%20b0c718fdb656445abd059bcb944703a4/Screenshot%209.png)

## ****HIVE scRNAseq Solution****

The HIVE collector is composed of a membrane with wells. The technology is based on pore properties and gravity sorting (more gentle approach) by the usage of vortex. It is possible to integrate sample storage in the workflow with freezing. The main advantage is the lack of specialized instrumentation, but we have a huge limit on the number of analyzable cells.