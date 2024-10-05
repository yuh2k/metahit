"""
detect the existence of depth.txt file before the echo "Coverage estimation completed successfully!"
remove all other files apart from depth.txt in the folder
"""
#!/usr/bin/env bash

#############################################################################################################
# Coverage Estimation Script using BBMap and jgi_summarize_bam_contig_depths.
#
# This script aligns WGS reads to the final assembly, sorts the BAM file by coordinate, and estimates the coverage.
#############################################################################################################

set -e  # Exit immediately if a command exits with a non-zero status
set -o pipefail  # Ensure that pipeline errors are not masked

# Display usage if not enough arguments
if [ "$#" -ne 8 ]; then
    echo "Usage: $0 -1 WGS_1.fastq.gz -2 WGS_2.fastq.gz -r final.contigs.fa -o output_dir"
    exit 1
fi

# Read in the arguments
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
    echo "Usage: $0 -1 WGS_1.fastq.gz -2 WGS_2.fastq.gz -r final.contigs.fa -o output_dir"
    exit 1
fi

mkdir -p "$output_dir"

# Increase Java heap memory to 4 GB (adjust as needed)
JAVA_HEAP="-Xmx4g"

# Run BBMap with increased memory
bbmap.sh in1="$reads_1" in2="$reads_2" ref="$ref" out="$output_dir/SG_map.sam" bamscript="$output_dir/bs.sh" $JAVA_HEAP

# Execute the generated BAM script
sh "$output_dir/bs.sh"

# Check if SAM file was created
if [[ ! -f "$output_dir/SG_map.sam" ]]; then
    echo "Error: SG_map.sam was not created."
    exit 1
fi

# Convert SAM to BAM, sort by coordinate, and index
samtools view -F 0x904 -bS "$output_dir/SG_map.sam" > "$output_dir/SG_map_UNSORTED.bam"
samtools sort "$output_dir/SG_map_UNSORTED.bam" -o "$output_dir/SG_map_sorted.bam"
samtools index "$output_dir/SG_map_sorted.bam"

# Run coverage estimation
./bin/metahit-scripts/jgi_summarize_bam_contig_depths --outputDepth "$output_dir/coverage.txt" --pairedContigs "$output_dir/pair.txt" "$output_dir/SG_map_sorted.bam"

echo "Coverage estimation completed successfully!"
