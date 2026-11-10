#!/bin/bash

#SBATCH --job-name=scRNA_index_QC_mapping_EAA
#SBATCH --time=48:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eacquah@uh.edu

set -euo pipefail

# Load required modules
module add java/25.0.4

# Software paths
FASTQC="/project/stuckert/eldrick/lab_6/FastQC/fastqc"
STAR="/project/stuckert/eldrick/lab_7/STAR-2.7.11b/bin/Linux_x86_64/STAR"

# Project directories
BASE="/project/stuckert/eldrick/Final-project"
DATA="$BASE/data"
REF="$BASE/references"
RESULTS="$BASE/results"

INDEX="$REF/STAR_index"
WHITELIST="$REF/3M-february-2018.txt"
INDEX_DONE="$INDEX/.index_complete"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "$RESULTS"


# STEP 1: Generate STAR genome index once

echo "Checking STAR genome index..."

if [[ ! -f "$INDEX_DONE" ]]; then

    echo "Generating STAR genome index..."

    # Remove incomplete index from previous runs.
    # Only do this when no other job uses it.
    if [[ -d "$INDEX" ]]; then
        rm -rf "$INDEX"
    fi

    mkdir -p "$INDEX"

    "$STAR" \
        --runThreadN "$THREADS" \
        --runMode genomeGenerate \
        --genomeDir "$INDEX" \
        --genomeFastaFiles "$REF/GRCz12ab.fa" \
        --sjdbGTFfile "$REF/GRCz12ab.gtf" \
        --sjdbOverhang 149 \
        --limitGenomeGenerateRAM 90000000000

    # Verify index files
    for FILE in genomeParameters.txt Genome SA SAindex; do
        if [[ ! -s "$INDEX/$FILE" ]]; then
            echo "ERROR: Missing index file: $FILE"
            exit 1
        fi
    done

    touch "$INDEX_DONE"

    echo "STAR genome indexing completed."

else
    echo "STAR genome index already exists."
fi


# STEP 2: Define and process samples sequentially

SAMPLES=("SRR24874370" "SRR24874372" "SRR24874373")
SAMPLE_GROUPS=("WT" "P23H" "P23H")

for i in "${!SAMPLES[@]}"; do

    SRR="${SAMPLES[$i]}"
    GROUP="${SAMPLE_GROUPS[$i]}"

    INPUT="$DATA/$GROUP"
    OUT="$RESULTS/$GROUP/$SRR"

    mkdir -p "$OUT"/{fastqc,STARsolo}

    BARCODE="$INPUT/${SRR}_2.fastq"
    CDNA="$INPUT/${SRR}_3.fastq"

    echo "Processing $SRR ($GROUP)"

    # Check input files
    for FILE in "$BARCODE" "$CDNA" "$WHITELIST"; do
        if [[ ! -s "$FILE" ]]; then
            echo "ERROR: Missing or empty file: $FILE"
            exit 1
        fi
    done

    # Check read counts
    N_CB=$(awk 'END {print NR/4}' "$BARCODE")
    N_CDNA=$(awk 'END {print NR/4}' "$CDNA")

    if [[ "$N_CB" != "$N_CDNA" ]]; then
        echo "ERROR: FASTQ read counts differ."
        exit 1
    fi

    echo "Reads per file: $N_CB"


# STEP 3: FastQC

echo "Running FastQC..."

if [[ ! -s "$OUT/fastqc/${SRR}_2_fastqc.zip" ||
      ! -s "$OUT/fastqc/${SRR}_3_fastqc.zip" ]]; then
   "$FASTQC" --threads "$THREADS" --outdir "$OUT/fastqc" "$BARCODE" "$CDNA"
else
   echo "FastQC already completed. Skipping."
fi


# STEP 4: STARsolo alignment

echo "Running STARsolo for $SRR..."

"$STAR" \
   --runThreadN "$THREADS" \
   --genomeDir "$INDEX" \
   --readFilesIn "$CDNA" "$BARCODE" \
   --outFileNamePrefix "$OUT/STARsolo/" \
   --outSAMtype None \
   --soloType CB_UMI_Simple \
   --soloCBstart 1 \
   --soloCBlen 16 \
   --soloUMIstart 17 \
   --soloUMIlen 12 \
   --soloBarcodeReadLength 0 \
   --soloCBwhitelist "$WHITELIST" \
   --soloFeatures Gene \
   --soloCellFilter EmptyDrops_CR \
   --soloUMIdedup 1MM_CR \
   --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts


# STEP 5: Verify count matrices

MATRIX="$OUT/STARsolo/Solo.out/Gene/filtered"

for FILE in matrix.mtx features.tsv barcodes.tsv; do
    if [[ ! -s "$MATRIX/$FILE" ]]; then
       echo "ERROR: Missing count matrix file: $FILE"
       exit 1
    fi
done

echo "Completed $SRR ($GROUP)"

done


echo "ALL THREE SAMPLES COMPLETED"
echo "STARsolo count matrices generated."
echo "Ready for downstream Seurat analysis."
