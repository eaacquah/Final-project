#!/bin/bash
set -euo pipefail   # safer bash



echo "=== QC + TRIM + MAP SCRIPT STARTED ==="


# Define directories

BASE=/project/stuckert/eldrick/Final-project
RAW=${BASE}/data/raw
TRIM=${BASE}/results/trimmed
QC=${BASE}/results/fastqc
MULTIQC=${BASE}/results/multiqc
ALIGN=${BASE}/results/alignments
REF=${BASE}/references
LOG=${BASE}/logs

mkdir -p $TRIM $QC $MULTIQC $ALIGN $LOG


# Load required modules

module load fastqc
module load multiqc
module load fastp
module load star


# 1. FASTQC on raw reads

echo "Running FastQC..."
fastqc $RAW/*.fastq* -o $QC


# 2. MultiQC summary report

echo "Running MultiQC..."
multiqc $QC -o $MULTIQC


# 3. Trim reads using fastp

echo "Running fastp trimming..."

for R1 in $RAW/*_1.fastq*
do
    R2=${R1/_1.fastq/_2.fastq}     # paired-end mate
    SAMPLE=$(basename $R1 | sed 's/_1.fastq//')

    fastp \
        -i $R1 \
        -I $R2 \
        -o $TRIM/${SAMPLE}_R1.trim.fastq.gz \
        -O $TRIM/${SAMPLE}_R2.trim.fastq.gz \
        --thread 8
done


# 4. Build STAR genome index (only once)

STAR_INDEX=${REF}/STAR_index
mkdir -p $STAR_INDEX

if [ ! -f "${STAR_INDEX}/SA" ]; then
    echo "Building STAR index..."
    STAR --runThreadN 8 \
         --runMode genomeGenerate \
         --genomeDir $STAR_INDEX \
         --genomeFastaFiles $REF/GRCz12.fa \
         --sjdbGTFfile $REF/GRCz12.gtf \
         --sjdbOverhang 99   # read_length - 1
fi


# 5. Align trimmed reads with STAR

echo "Running STAR alignment..."

for R1 in $TRIM/*_R1.trim.fastq.gz
do
    R2=${R1/_R1.trim.fastq.gz/_R2.trim.fastq.gz}
    SAMPLE=$(basename $R1 | sed 's/_R1.trim.fastq.gz//')

    STAR --runThreadN 8 \
         --genomeDir $STAR_INDEX \
         --readFilesIn $R1 $R2 \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts \
         --outFileNamePrefix $ALIGN/${SAMPLE}_
done

echo "=== QC + TRIM + MAP SCRIPT FINISHED ==="
