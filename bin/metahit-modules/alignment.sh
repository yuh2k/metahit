#!/usr/bin/env bash

#############################################################################################################
# Alignment script using BWA MEM and SamTools.
#
# Aligns reads to the final assembly and processes the output using SamTools.
#############################################################################################################

# Usage
echo "Usage: metahit alignment -f final_assembly.fasta -1 HIC_1_dedup.fastq.gz -2 HIC_2_dedup.fastq.gz -o output_dir"

# Parameters
fasta="final_assembly.fasta"
reads_1="HIC_1_dedup.fastq.gz"
reads_2="HIC_2_dedup.fastq.gz"阿
output_dir=".output/final/"

# Run BWA MEM to align reads
bwa mem -5SP $fasta $reads_1 $reads_2 > $output_dir/MAP.sam

# Convert SAM to BAM, filtering out unwanted flags
samtools view -F 0x904 -bS $output_dir/MAP.sam > $output_dir/MAP_UNSORTED.bam

# Sort BAM file by read names
samtools sort -n $output_dir/MAP_UNSORTED.bam -o $output_dir/MAP_SORTED.bam

# Clean up
rm $output_dir/MAP.sam $output_dir/MAP_UNSORTED.bam

echo "Alignment completed successfully!"