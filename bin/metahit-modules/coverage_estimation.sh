#!/usr/bin/env bash

# Coverage Estimation Script using BBMap and jgi_summarize_bam_contig_depths.

set -e
set -o pipefail

# Function to print informational messages
function echo_info() {
    echo -e "\033[1;34m[INFO]\033[0m $1"
}

# Function to print error messages
function echo_error() {
    echo -e "\033[1;31m[ERROR]\033[0m $1" >&2
}

# Display usage if not enough arguments
if [ "$#" -lt 8 ]; then
    echo "Usage: $0 -1 WGS_1.fastq.gz -2 WGS_2.fastq.gz -r final.contigs.fa -o output_dir [-m Java_heap_memory]"
    exit 1
fi

# Set default Java heap memory to 80% of available memory
default_memory=$(( $(free -m | awk '/^Mem:/{print $2}') * 80 / 100 ))
JAVA_HEAP="-Xmx${default_memory}m"

# Read in the arguments
while getopts "p:1:2:r:o:m:" opt; do
    case $opt in
        p) path=$OPTARG ;;
        1) reads_1=$OPTARG ;;
        2) reads_2=$OPTARG ;;
        r) ref=$OPTARG ;;
        o) output_dir=$OPTARG ;;
        m) JAVA_HEAP="-Xmx$OPTARG" ;;
        *) echo "Invalid option"; exit 1 ;;
    esac
done

# Check if required parameters are set
if [[ -z "$reads_1" || -z "$reads_2" || -z "$ref" || -z "$output_dir" ]]; then
    echo_error "Missing required parameters"
    echo "Usage: $0 -1 WGS_1.fastq.gz -2 WGS_2.fastq.gz -r final.contigs.fa -o output_dir [-m Java_heap_memory]"
    exit 1
fi


# Run BBMap with memory setting and capture output for debugging
echo_info "Running BBMap..."
${path}/external/bbmap/bbmap.sh in1="$reads_1" in2="$reads_2" ref="$ref" out="$output_dir/SG_map.sam" bamscript="$output_dir/bs.sh" $JAVA_HEAP > "$output_dir/bbmap.log" 2>&1

# Check if BBMap generated the SAM file successfully
if [[ ! -f "$output_dir/SG_map.sam" ]]; then
    echo_error "BBMap failed to generate SG_map.sam. Check $output_dir/bbmap.log for details."
    exit 1
fi

# Run jgi_summarize_bam_contig_depths and log output
echo_info "Running jgi_summarize_bam_contig_depths..."
${path}/bin/metahit-scripts/jgi_summarize_bam_contig_depths \
    --outputDepth "$output_dir/coverage.txt" \
    --pairedContigs pair.txt \
    "$output_dir/SG_map.sam" > "$output_dir/jgi_summarize.log" 2>&1

# Check if coverage.txt was created
if [[ -f "$output_dir/coverage.txt" ]]; then
    echo_info "Coverage.txt generated successfully at $output_dir/estimation/coverage.txt"
else
    echo_error "Failed to generate coverage.txt. Check $output_dir/jgi_summarize.log for details."
    exit 1
fi
find "$output_dir" -type f ! -name "coverage.txt" -exec rm -f {} +

