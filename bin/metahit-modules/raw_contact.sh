#!/usr/bin/env bash
echo "[FREE MEMORY]: $free_mem"
# Check if the correct number of arguments is passed
if [ "$#" -lt 6 ]; then
    echo "Usage: $0 --bam <BAM file> --fasta <FASTA file> --out <output directory> --enzyme <enzyme name> [optional parameters]"
    exit 1
fi

# Parse input arguments
BAM=""
FASTA=""
OUTDIR=""
ENZYME=""
METACC_MIN_SIGNAL=1
METACC_MIN_LEN=1000
METACC_MIN_MAPQ=30
METACC_MIN_MATCH=30

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -p) path=$2; shift 2;;
        --bam) BAM="$2"; shift 2;;
        --fasta) FASTA="$2"; shift 2;;
        --out) OUTDIR="$2"; shift 2;;
        --enzyme) ENZYME="$2"; shift 2;;
        --metacc-min-signal) METACC_MIN_SIGNAL="$2"; shift 2;;
        --metacc-min-len) METACC_MIN_LEN="$2"; shift 2;;
        --metacc-min-mapq) METACC_MIN_MAPQ="$2"; shift 2;;
        --metacc-min-match) METACC_MIN_MATCH="$2"; shift 2;;
        *) echo "[ERROR] Unknown parameter: $1"; exit 1;;
    esac
done

# Check if required arguments are provided
if [ -z "$BAM" ] || [ -z "$FASTA" ] || [ -z "$OUTDIR" ] || [ -z "$ENZYME" ]; then
    echo "[ERROR] Missing required arguments."
    echo "Usage: $0 --bam <BAM file> --fasta <FASTA file> --out <output directory> --enzyme <enzyme name> [optional parameters]"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"


# Run the Python script to generate raw contact matrices with the provided parameters
python ${path}/bin/metahit-scripts/raw_contact.py \
    --BAM "$BAM" \
    --FASTA "$FASTA" \
    --OUTDIR "$OUTDIR" \
    --enzyme "$ENZYME" \
    --metacc-min-signal "$METACC_MIN_SIGNAL" \
    --metacc-min-len "$METACC_MIN_LEN" \
    --metacc-min-mapq "$METACC_MIN_MAPQ" \
    --metacc-min-match "$METACC_MIN_MATCH"

# Check if the Python script executed successfully
if [ $? -eq 0 ]; then
    echo "[INFO] Raw contact matrix generation completed successfully."
else
    echo "[ERROR] Raw contact matrix generation failed."
    exit 1
fi
