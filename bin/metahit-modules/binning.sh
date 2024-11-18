#!/usr/bin/env bash

# Activate your environment if necessary
# source activate metahit_env

# Define input parameters
OUTPUT_DIR="output/bins"
FASTA_FILE="/path/to/final_assembly.fasta"  # Replace with the actual path
BAM_FILE="/path/to/sorted_map.bam"          # Replace with the actual path
COVERAGE_FILE="/path/to/coverage.txt"       # Replace with the actual path
ENZYMES="HindIII"                           # Replace with your enzymes, e.g., "HindIII MboI"
NUM_GENE=100                                # Replace with the actual number or leave empty
SEED=42                                     # Replace with desired seed or leave empty

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run the bin refinement script
echo "[INFO] Running bin refinement..."
python bin_refinement.py \
  --bam "$BAM_FILE" \
  --enzyme $ENZYMES \
  --fasta "$FASTA_FILE" \
  --coverage "$COVERAGE_FILE" \
  --output "$OUTPUT_DIR" \
  ${NUM_GENE:+--num_gene "$NUM_GENE"} \
  ${SEED:+--seed "$SEED"}

echo "[INFO] Binning processes completed."
