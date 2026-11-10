# Final Project Workflow
Author: Eldrick Acquah 
Project: Transcriptomic changes in the P23H Zebrafish model of Retinitis Pigmentosa 


## Overview
This workflow describes the complete computational pipeline used to process, analyze, and interpret single‑cell RNA‑seq data from WT and P23H zebrafish samples. The workflow consists of four main scripts:

1. install_and_download.sh 
   Downloads all required datasets, reference genome files, and raw FASTQ reads.

2. cell_barcode_whitelist_matching.sh
   Performs sequence read structure and barcode whitelist matching compatibility checks

3. qc_index_map.sh
   Performs quality control, read trimming, genome indexing, and alignment.

4. analysis.R 
   Conducts differential gene expression analysis using DESeq2 (pseudo‑bulk) or Seurat (single‑cell).

Scripts 1 to 3 will run on the Sabine HPC cluster using SLURM job submission, whiles script 3 (analysis.R) will run on RStudio


# Repository Structure
```
Final-project/
├── README.md
├── .gitignore
├── Workflow.md
├── config/
│   └── samples.tsv
├── scripts/
│   ├── install_and_download.sh
│   ├── qc_trim_map.sh
│   └── analysis.R
├── references/
│   ├── GRCz12.fa
│   └── GRCz12.gtf
├── data/
│   └── raw/               # FASTQ files (ignored by Git)
├── results/
│   ├── fastqc/
│   ├── multiqc/
│   ├── trimmed/
│   ├── alignments/
│   └── counts/
└── logs/
```

Large files (FASTQ, BAM, STAR index, MultiQC reports) are excluded from version control using `.gitignore`.


## Pipeline Summary

# 1. Data Acquisition & Setup
Script: "scripts/install_and_download.sh"

This script:
- Loads required modules on Sabine  
- Downloads the GRCz12 zebrafish genome and GTF annotation from NCBI/Ensembl
- Downloads raw FASTQ files from GEO/SRA using "prefetch" + "fasterq-dump"  
- Stores all raw data in "data/raw/"  
- Logs all actions to "logs/install_manifest.txt"

This step ensures the project has all necessary inputs before processing.



# 2. Quality Control, Trimming, and Alignment
Script: "scripts/q_trim_map.sh"

This script performs:

# FastQC
Runs FastQC on all raw FASTQ files and stores reports in "results/fastqc/".

# MultiQC
Aggregates FastQC results into a single HTML report in "results/multiqc/".

# STAR Genome Indexing
Builds a STAR index using:
- "references/GRCz12.fa"
- "references/GRCz12.gtf"

Index stored in "references/STAR_index/".

# STAR Alignment
Aligns trimmed reads to the GRCz12 genome using STAR:
- Produces sorted BAM files 
- Generates gene‑level counts via "--quantMode GeneCounts" 
- Stores alignment outputs in "results/alignments/"

This step produces the count tables needed for downstream differential expression.



# 3. Differential Expression Analysis
Script: "scripts/analysis.R

Seurat-based  workflow for cell‑level analysis:
1. Load cell‑by‑gene matrices (Cell Ranger or STARsolo) 
2. Perform QC filtering 
3. Normalize (SCTransform or LogNormalize) 
4. Dimensionality reduction (PCA, UMAP) 
5. Clustering 
6. Marker gene identification 
7. Save Seurat object to "results/seurat_object.rds"

This mode is appropriate for true single‑cell datasets.


# SLURM Usage
The first three scripts will be submitted to Sabine using:

```bash
sbatch scripts/install_and_download.sh
sbatch scripts/qc_trim_map.sh
```

The following cab be adjusted:

- `--cpus-per-task`
- `--mem`
- `--time`
- `--nodes`

based on dataset size.


# Reproducibility
This workflow ensures reproducibility through:

- GitHub version control 
- `.gitignore` preventing large data from being committed 
- Modular scripts 
- Logged outputs in `logs/` 



# Notes
- Raw FASTQ and BAM files are **never** committed to GitHub.  
- STAR index files are also excluded from version control.  
- All scripts are designed to be portable and cluster‑friendly.  
- The workflow supports both bulk‑style and single‑cell analysis paths.



## Detailed scRNA-Seq Processing Workflow
1. Data Acquisition and Subsampling
Three SRA accessions were selected for the comparison of WT and P23H transgenic zebrafish retinas:
- WT: SRR24874370
- P23H: SRR24874372
- P23H: SRR24874373
SRA Toolkit was used to obtain FASTQ reads. A subset of 10 million sequencing spots per accession was prepared and used. 
The FASTQ files were stored under data/WT/ and data/P23H/.

2. FASTQ Structure and Barcode Validation
The _2.fastq files were assigned as barcode/UMI reads.
The _3.fastq files were assigned as cDNA reads.
The expected 10x Chromium 3' v3 layout was used:
- Cell barcode: bases 1–16
- UMI: bases 17–28
- Whitelist: 3M-february-2018.txt
Barcode-whitelist matching was checked on a read subset. The read pairing and counts were also checked before mapping.

3. Read Quality Assessment
FastQC was run separately on barcode/UMI and cDNA FASTQ files. Quality reports were directed to each sample's fastqc/ folder.
Barcode/UMI positions were preserved for STARsolo processing.

4. Genome Reference and Indexing
The following files were present in references/:
- GRCz12ab.fa
- GRCz12ab.gtf
- 3M-february-2018.txt
STAR v2.7.11b was configured to generate a genome index using the GRCz12ab FASTA and GTF files. The indexing command used
eight threads and --sjdbOverhang 149.
Index generation occured once for all samples.

5. STARsolo Alignment and UMI Counting
STARsolo was configured with these settings:
--soloType CB_UMI_Simple
--soloCBstart 1
--soloCBlen 16
--soloUMIstart 17
--soloUMIlen 12
--soloBarcodeReadLength 0
--soloCBwhitelist references/3M-february-2018.txt
--soloFeatures Gene
--soloCellFilter EmptyDrops_CR
--soloUMIdedup 1MM_CR
--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
--outSAMtype None

The cDNA FASTQ was supplied before the barcode/UMI FASTQ in --readFilesIn. The --outSAMtype None option disabled BAM
output. STARsolo produced UMI count matrices.

Output for each sample:
STARsolo/Solo.out/Gene/
├── raw/
│   ├── matrix.mtx
│   ├── features.tsv
│   └── barcodes.tsv
└── filtered/
    ├── matrix.mtx
    ├── features.tsv
    └── barcodes.tsv

The filtered/ directory was the input for Seurat.

6. Downstream Seurat Analysis
Seurat was selected for downstream processing in R. The analysis included cell-level QC, normalization, dimensionality reduction, clustering, retinal cell-type annotation, and WT-versus-P23H expression comparisons within cell populations.
Microglia, rods, cones, and Müller glia were priority cell types.
