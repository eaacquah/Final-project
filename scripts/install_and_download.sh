#!/bin/bash

#SBATCH --job-name=Install_and_download_script_EAA
#SBATCH --mem=55G
#SBATCH -N 1 -n 8
#SBATCH -t 5:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eaacquah@uh.edu 


echo "=== INSTALL & DOWNLOAD SCRIPT STARTED ==="


# Setting project directories

BASE=/project/stuckert/eldrick/Final-project
REF=${BASE}/references        # reference genome + annotation
RAW=${BASE}/data             # raw FASTQ files
LOG=${BASE}/logs              # logs for reproducibility

mkdir -p $REF $RAW $LOG       # create directories if missing



# Loading required modules

module load sratoolkit/3.4.1         # for prefetch + fasterq-dump
module load python/3.12               # required by some tools

echo "Modules loaded." >> ${LOG}/install_manifest.txt


# 1. Download GRCz12 genome + annotation (FASTA + GTF)

cd "$REF"

echo "Downloading GRCz12ab genome and annotation..."

# Download genome FASTA
wget -c -O GRCz12ab.fa.gz \
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/052/040/795/GCF_052040795.1_GRCz12ab/GCF_052040795.1_GRCz12ab_genomic.fna.gz

# Download annotation GTF
wget -c -O GRCz12ab.gtf.gz \
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/052/040/795/GCF_052040795.1_GRCz12ab/GCF_052040795.1_GRCz12ab_genomic.gtf.gz

# Unzip
gunzip GRCz12ab.fa.gz
gunzip GRCz12ab.gtf.gz

echo "GRCz12ab genome + GTF downloaded." \
>> "${LOG}/install_manifest.txt"


# 2. Download FASTQ files from GEO/SRA

cd $RAW

echo "Downloading FASTQ files..."

# Create sample directories
mkdir -p "$RAW/WT"
mkdir -p "$RAW/P23H"

# WT sample
WT_ACCESSIONS=("SRR24874370")

# P23H samples
P23H_ACCESSIONS=("SRR24874372" "SRR24874373")

# Download WT reads
echo "Downloading WT samples..."

for SRR in "${WT_ACCESSIONS[@]}"; do
    echo "Processing $SRR..."

    fastq-dump --split-files --skip-technical -X 10000000 "$SRR" --outdir "$RAW/WT"

done

# Download P23H reads
echo "Downloading P23H samples..."

for SRR in "${P23H_ACCESSIONS[@]}"; do
    echo "Processing $SRR..."

    fastq-dump \
        --split-files \
        --skip-technical \
        -X 10000000 \
        "$SRR" \
        --outdir "$RAW/P23H"

done

echo "All downloads completed!" 

echo "=== INSTALL & DOWNLOAD SCRIPT FINISHED ==="

