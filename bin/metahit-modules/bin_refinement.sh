#!/usr/bin/env bash
free_mem=$(free -h | awk '/^Mem:/ {print $4}')
echo "[FREE MEMORY]: $free_mem"
# Script to perform bin refinement for Metahit

# Display usage if not enough arguments
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <FASTA> <BAM> <OUTDIR> [options]"
    echo "Required arguments:"
    echo "  <FASTA>          Reference fasta sequence"
    echo "  <BAM>            Input bam file in query order"
    echo "  <OUTDIR>         Output directory"
    echo "Optional arguments:"
    echo "  -e, --enzyme        Case-sensitive enzyme name. Use multiple times for multiple enzymes"
    echo "  --metacc-min-len    Minimum acceptable contig length [1000]"
    echo "  --metacc-min-signal Minimum acceptable Hi-C signal [2]"
    echo "  --metacc-min-mapq   Minimum acceptable mapping quality [30]"
    echo "  --threads          The number of threads. Default is 30."
    echo "  --seed             Random seed"
    exit 1
fi

FASTA=$1
BAM=$2
OUTDIR=$3
path=$4
shift 4

eval "$(conda shell.bash hook)"
conda activate checkm2
# Path to the bin refinement Python script
BIN_REFINEMENT_SCRIPT="${path}/bin/metahit-scripts/bin_refinement.py"

# Create output directory if it does not exist
if [ ! -d "$OUTDIR" ]; then
    mkdir -p "$OUTDIR"
fi

# Run the bin refinement script with provided arguments
python "$BIN_REFINEMENT_SCRIPT" --FASTA "$FASTA" --BAM "$BAM" --OUTDIR "$OUTDIR" "$@"

# Check if the bin refinement process succeeded
if [ $? -ne 0 ]; then
    echo "Error: Bin refinement process failed."
    exit 1
fi

# Inform that the process has been completed successfully
echo "Bin refinement completed successfully."
conda deactivate