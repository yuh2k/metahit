#!/usr/bin/env bash

# Usage: ./normalization.sh <command> [options]
# Example: ./normalization.sh normcc -c contig_info.csv -m contact_matrix.npz -o output --refinement_method metacc

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <command> [options]"
    echo "Available normalization commands: raw, normcc, hiczin, bin3c, metator, fastnorm"
    echo "Optional refinement method: --refinement_method <metacc|bin3c|imputecc>"
    exit 1
fi

COMMAND=$1
shift

REFINEMENT_METHOD=""
OUTDIR=""

# Parse options
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --refinement_method)
            REFINEMENT_METHOD="$2"
            shift 2
            ;;
        -o|--output)
            OUTDIR="$2"
            shift 2
            ;;
        *)
            NORMALIZATION_ARGS+=("$1")
            shift
            ;;
    esac
done

# Run the normalization command
./bin/metahit-scripts/normalization.py "$COMMAND" --contig_file "$CONTIG_FILE" --contact_matrix_file "$CONTACT_MATRIX_FILE" --output "$OUTDIR" "${NORMALIZATION_ARGS[@]}"


# Check if refinement is requested
if [ -n "$REFINEMENT_METHOD" ]; then
    echo "[INFO] Running refinement using method: $REFINEMENT_METHOD"
    
    CONTIG_FILE="${OUTDIR}/contig_info.csv"
    HIC_MATRIX="${OUTDIR}/denoised_contact_matrix_${COMMAND}.npz"
    
    if [ ! -f "$CONTIG_FILE" ] || [ ! -f "$HIC_MATRIX" ]; then
        echo "[ERROR] Required files for refinement not found."
        exit 1
    fi

    case $REFINEMENT_METHOD in
        metacc)
            ./bin/metahit-scripts/bin_refinement.py --method metacc --contig_file "$CONTIG_FILE" --hic_matrix "$HIC_MATRIX" --outdir "$OUTDIR"
            ;;
        bin3c)
            ./bin/metahit-scripts/bin_refinement.py --method bin3c --contig_file "$CONTIG_FILE" --hic_matrix "$HIC_MATRIX" --outdir "$OUTDIR" --fasta "${OUTDIR}/reference.fasta"
            ;;
        imputecc)
            ./bin/metahit-scripts/bin_refinement.py --method imputecc --contig_file "$CONTIG_FILE" --hic_matrix "$HIC_MATRIX" --outdir "$OUTDIR"
            ;;
        *)
            echo "[ERROR] Invalid refinement method: $REFINEMENT_METHOD"
            exit 1
            ;;
    esac
fi
