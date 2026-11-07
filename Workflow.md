# Final Project Workflow
Author: Eldrick Acquah  
Project: Transcriptomic changes in the P23H Zebrafish model of Retinitis Pigmentosa 


# Overview
This workflow describes the complete computational pipeline used to process, analyze, and interpret single‑cell RNA‑seq data from WT and P23H zebrafish samples. The workflow is fully reproducible, modular, and tracked through GitHub. It consists of three main scripts:

1. install_and_download.sh  
   Downloads all required datasets, reference genome files, and raw FASTQ reads.

2. qc_trim_map.sh
   Performs quality control, read trimming, genome indexing, and alignment.

3. analysis.R  
   Conducts differential gene expression analysis using DESeq2 (pseudo‑bulk) or Seurat (single‑cell).

Scripts 1 and 2 will run on the Sabine HPC cluster using SLURM job submission, whiles script 3 (analysis.R) will run on RStudio.


# Repository Structure
```
Final-project/
├── README.md
├── .gitignore
├── workflow.md
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

# fastp Trimming
Performs adapter trimming and quality filtering for paired‑end reads.  
Outputs trimmed reads to "results/trimmed/".

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
Script: "scripts/analysis.R"

Two analysis modes are supported:

# A. DESeq2 (Pseudo‑Bulk)
Uses STAR’s `ReadsPerGene.out.tab` files.

Workflow:

1. Load count tables 
2. Build a DESeq2 dataset 
3. Model WT vs P23H differences 
4. Generate differential expression results 
5. Save results to "results/DESeq2_results_P23H_vs_WT.csv"

This mode is appropriate when aggregating counts per sample.

---

# B. Seurat (Single‑Cell)
Optional workflow for cell‑level analysis:

1. Load cell‑by‑gene matrices (Cell Ranger or STARsolo) 
2. Perform QC filtering 
3. Normalize (SCTransform or LogNormalize) 
4. Dimensionality reduction (PCA, UMAP) 
5. Clustering 
6. Marker gene identification 
7. Save Seurat object to "results/seurat_object.rds"

This mode is appropriate for true single‑cell datasets.


# SLURM Usage
The first two scripts will be submitted to Sabine using:

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

