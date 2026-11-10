#!/bin/bash 

# Download reference for matching 
wget -O 3M-february-2018.txt.gz \
https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/3M-february-2018.txt.gz

# Verify it is a valid gzip file
gzip -t 3M-february-2018.txt.gz

# Decompress
gunzip 3M-february-2018.txt.gz

# Count barcodes
wc -l 3M-february-2018.txt

# Inspect the first five barcodes
head -n 5 3M-february-2018.txt

cd /project/stuckert/eldrick/Final-project


# Performing whitelist matching with python script

python - <<'PY'
from pathlib import Path

# 10x Genomics v3 barcode whitelist
whitelist_path = Path("references/3M-february-2018.txt")

with open(whitelist_path) as f:
    whitelist = {line.strip() for line in f}

# Sample information
samples = {
    "WT": ["SRR24874370"],
    "P23H": ["SRR24874372", "SRR24874373"]
}

max_reads = 100000

print(f"{'Genotype':<12} {'Accession':<16} {'FASTQ':<8} "
      f"{'Reads':<12} {'Matches':<12} {'Match (%)':<12}")
print("-" * 75)

for genotype, accessions in samples.items():
    for srr in accessions:
        for read_number in (2, 3):

            fastq = Path(f"data/{genotype}/{srr}_{read_number}.fastq")

            if not fastq.exists():
                print(f"File not found: {fastq}")
                continue

            total = 0
            matches = 0

            with open(fastq) as f:
                while total < max_reads:
                    header = f.readline()
                    if not header:
                        break

                    sequence = f.readline().strip()
                    f.readline()  # Plus line
                    f.readline()  # Quality scores

                    barcode = sequence[:16]

                    if barcode in whitelist:
                        matches += 1

                    total += 1

            percentage = (matches / total * 100) if total else 0

            print(f"{genotype:<12} {srr:<16} "
                  f"R{read_number:<7} {total:<12} "
                  f"{matches:<12} {percentage:<12.2f}")
PY

