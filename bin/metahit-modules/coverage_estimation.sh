
#!/usr/bin/env bash
#############################################################################################################
# Coverage Estimation Script using BBMap and jgi_summarize_bam_contig_depths.
#
# This script aligns WGS reads to the final assembly, sorts the BAM file by read name for normalization,
# then sorts by coordinate for indexing, and estimates the coverage.
#############################################################################################################

set -e  # Exit immediately if a command exits with a non-zero status
set -o pipefail  # Ensure that pipeline errors are not masked

# Function to print informational messages
function echo_info() {
    echo -e "\033[1;34m[INFO]\033[0m $1"
}

# Function to print error messages
function echo_error() {
    echo -e "\033[1;31m[ERROR]\033[0m $1" >&2
}

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
    echo_error "Missing required parameters"
    echo "Usage: $0 -1 WGS_1.fastq.gz -2 WGS_2.fastq.gz -r final.contigs.fa -o output_dir"
    exit 1
fi

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Increase Java heap memory to 4 GB (adjust as needed)
JAVA_HEAP="-Xmx4g"

# Run BBMap with increased memory
./external/bbmap/bbmap.sh in1="$reads_1" in2="$reads_2" ref="$ref" out="$output_dir/SG_map.sam" bamscript="$output_dir/bs.sh" $JAVA_HEAP

bin/metahit-scripts/jgi_summarize_bam_contig_depths --outputDepth output/estimation/coverage.txt --pairedContigs pair.txt $output_dir/SG_map.sam

find "$output_dir" -type f ! -name "coverage.txt" -exec rm -f {} +


# Check if coverage.txt was created
if [[ -f "$output_dir/coverage.txt" ]]; then
    echo_info "Coverage.txt generated successfully at $output_dir/coverage.txt"
else
    echo_error "Failed to generate coverage.txt. Check $output_dir/jgi_summarize.log for details."
    exit 1
fi

# Final success message
echo_info "Coverage estimation completed successfully."
