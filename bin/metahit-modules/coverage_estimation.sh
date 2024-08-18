#!/usr/bin/env bash

#############################################################################################################
# Coverage Estimation Script using BBMap and jgi_summarize_bam_contig_depths.
#
# This script aligns WGS reads to the final assembly, sorts the BAM file, and estimates the coverage.
#############################################################################################################

echo "Usage: metahit coverage -1 WGS_1.fastq.gz -2 WGS_2.fastq.gz -r final.contigs.fa -o output_dir"

# Input parameters
reads_1=""
reads_2=""
ref=""
output_dir=""

# Parse command line arguments
while getopts "1:2:r:o:" opt; do
  case $opt in
    1) reads_1=$OPTARG ;;
    2) reads_2=$OPTARG ;;
    r) ref=$OPTARG ;;
    o) output_dir=$OPTARG ;;
    *) echo "Invalid option"; exit 1 ;;
  esac
done

# Check if required parameters are set
if [[ -z "$reads_1" || -z "$reads_2" || -z "$ref" || -z "$output_dir" ]]; then
  echo "Missing required parameters"
  echo "Usage: metahit coverage -1 WGS_1.fastq.gz -2 WGS_2.fastq.gz -r final.contigs.fa -o output_dir"
  exit 1
fi

mkdir -p $output_dir
bbmap.sh in1=$reads_1 in2=$reads_2 ref=$ref out=$output_dir/SG_map.sam bamscript=$output_dir/bs.sh
sh $output_dir/bs.sh
./bin/metahit-scripts/jgi_summarize_bam_contig_depths --outputDepth $output_dir/coverage.txt --pairedContigs $output_dir/pair.txt $output_dir/SG_map_sorted.bam

echo "Coverage estimation completed successfully!"
