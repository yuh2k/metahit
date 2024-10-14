#!/usr/bin/env bash

set -e
set -o pipefail
set -x

# Define tool paths
BWA_PATH="./external/bin/bwa"
SAMTOOLS_PATH="./external/bin/samtools"

# Define input and output files
REFERENCE="output/assembly/final_assembly.fasta"
OUTPUT_DIR="output/alignment"
READS_1="output/readqc/hic/final_pure_reads_1.fastq"
READS_2="output/readqc/hic/final_pure_reads_2.fastq"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Index the reference genome
echo "Indexing reference genome with BWA..."
$BWA_PATH index "$REFERENCE"

# Align reads with BWA MEM
echo "Aligning reads with BWA MEM..."
$BWA_PATH mem -t 4 "$REFERENCE" "$READS_1" "$READS_2" > "$OUTPUT_DIR/map.sam"

# Convert SAM to BAM
echo "Converting SAM to BAM..."
$SAMTOOLS_PATH view -F 0x904 -bS "$OUTPUT_DIR/map.sam" > "$OUTPUT_DIR/unsorted_map.bam"

# Sort BAM by read name
echo "Sorting BAM by read name..."
$SAMTOOLS_PATH sort -n "$OUTPUT_DIR/unsorted_map.bam" -o "$OUTPUT_DIR/sorted_map.bam"


echo "Alignment completed successfully!"
