#!/bin/bash

# Help function to display usage
function show_help {
    echo "Usage: $0 --bam BAM_FILE --fasta FASTA_FILE --out OUTPUT_DIR"
    echo ""
    echo "Required arguments:"
    echo "  --bam         Path to the BAM file containing Hi-C reads"
    echo "  --fasta       Path to the FASTA file containing contig sequences"
    echo "  --out         Output directory to save Hi-C contact matrix (.npz) and contig info (.csv)"
    echo ""
    echo "Example usage:"
    echo "$0 --bam output/alignment/MAP_SORTED.bam --fasta output/assembly/final_assembly.fasta --out output/downstream"
}

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --bam) BAM_FILE="$2"; shift ;;
        --fasta) FASTA_FILE="$2"; shift ;;
        --out) OUTPUT_DIR="$2"; shift ;;
        --help) show_help; exit 0 ;;
        *) echo "Unknown parameter passed: $1"; show_help; exit 1 ;;
    esac
    shift
done

# Check if required arguments are provided
if [ -z "$BAM_FILE" ] || [ -z "$FASTA_FILE" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments."
    show_help
    exit 1
fi

# Print the received arguments
echo "BAM file: $BAM_FILE"
echo "FASTA file: $FASTA_FILE"
echo "Output directory: $OUTPUT_DIR"

# Run the downstream.py script with the provided arguments
python bin/metahit-scripts/raw_contig.py --bam "$BAM_FILE" --fasta "$FASTA_FILE" --out "$OUTPUT_DIR"
