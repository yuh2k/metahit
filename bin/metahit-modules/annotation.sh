#!/usr/bin/env bash

echo "Running: gtdbtk classify_wf"

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -p)
            path=$2
            shift 2
            ;;
        *)
            break
            ;;
    esac
done

# Check for required argument
if [ -z "$path" ]; then
    echo "Error: MetaHit path (-p) not provided."
    exit 1
fi

# Create output directory
mkdir -p "${path}/output/annotation"


export GTDBTK_DATA_PATH="/home/hh/Documents/gtdbtk_r220_data.tar.gz.part_aa"
if [ ! -d "$GTDBTK_DATA_PATH" ]; then
    echo "Error: GTDB-Tk data path is invalid or not set correctly. Please verify."
    exit 1
fi

# Activate the GTDB-Tk environment
echo "[INFO] Activating GTDB-Tk environment 'gtdbtk-2.4.0'..."
eval "$(conda shell.bash hook)"
conda activate gtdbtk-2.4.0
if [ $? -ne 0 ]; then
    echo "Error: Could not activate Conda environment 'gtdbtk-2.4.0'."
    exit 1
fi

# Run GTDB-Tk annotation
echo "[INFO] Running GTDB-Tk classify_wf..."
gtdbtk classify_wf --skip_ani_screen "$@" 
if [ $? -ne 0 ]; then
    echo "Error: Annotation (GTDB-Tk) step failed."
    exit 1
fi

echo "Annotation (GTDB-Tk) step completed successfully."
