# This is a workflow for my final project for Bioinformatics Fall 2026
## Single-Cell RNA-Seq Analysis of P23H Zebrafish Retinas

## Overview
This repository documents a bioinformatics workflow for processing single-cell RNA sequencing (scRNA-seq) data from wild-type (WT) and P23H transgenic zebrafish retinas. 
The analysis focused on retinal cell populations and transcriptional responses associated with photoreceptor degeneration, particularly microglial inflammation.

## Dataset
Three NCBI SRA accessions were selected:
Accession	Genotype
SRR24874370	WT
SRR24874372	P23H
SRR24874373	P23H

The reads were were obtained with the 10x Genomics Chromium 3' v3 platform and treated as such.
A subset of 10 million sequencing spots per accession was extracted and used for the analysis due to the large size of the 
complete dataset. 
The GRCz12ab zebrafish genome and matching GTF annotation were used as references.

## Workflow
1. SRA reads were downloaded and converted to FASTQ.
2. FASTQ read structure and  barcode-whitelist compatibility were inspected.
3. FastQC was used to assess sequence quality.
4. STAR genome indexing was configured for GRCz12ab.
5. STARsolo was used for sequence alignment and UMI counting.
6. Seurat was selected for downstream single-cell analysis.

See [Workflow.md](Workflow.md) for commands, expected outputs, and the status of each stage.

