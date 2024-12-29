#!/usr/bin/env bash

echo "[INFO] Running Virus-Host Interaction Analysis"

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --BIN) bin_file=$2; shift ;;
        --viral-contig) viral_contig_file=$2; shift ;;
        --contact) contact_file=$2; shift ;;
        --OUTDIR) out_dir=$2; shift ;;
        -p) metahit_path=$2; shift ;;
        -t) threads=$2; shift ;;
        -m) memory=$2; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Check for required arguments
if [ -z "$bin_file" ] || [ -z "$viral_contig_file" ] || [ -z "$contact_file" ] || [ -z "$out_dir" ]; then
    echo "[ERROR] Missing required arguments."
    echo "Usage: virus_host_interaction.sh --BIN <BIN_FILE> --viral-contig <VIRAL_CONTIG> --contact <CONTACT_MATRIX> --OUTDIR <OUTPUT_DIR> [-p <PATH>] [-t <THREADS>] [-m <MEMORY>]"
    exit 1
fi

# Create output directory
mkdir -p "$out_dir"

# Activate the metahit environment
eval "$(conda shell.bash hook)"
if ! conda activate metahit_env; then
    echo "[ERROR] Could not activate Conda environment 'metahit_env'."
    exit 1
fi

# Run the virus_host_interaction Python script
python "${metahit_path}/bin/metahit-modules/virus_host_interaction.py" \
    --BIN "$bin_file" \
    --viral-contig "$viral_contig_file" \
    --contact "$contact_file" \
    --OUTDIR "$out_dir" \
    -t "$threads" \
    -m "$memory"

if [ $? -ne 0 ]; then
    echo "[ERROR] Virus-Host Interaction step failed."
    exit 1
fi

echo "[INFO] Virus-Host Interaction completed successfully."
