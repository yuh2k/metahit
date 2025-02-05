#!/usr/bin/env bash
# genomad.sh - Pipeline for plasmid/virus sequence identification using genomad.
# This script automatically activates the genomad conda environment (assumed to be named "genomad_env")
# and allows the user to specify the --splits parameter (default is 8).

free_mem=$(free -h | awk '/^Mem:/ {print $4}')
echo "[FREE MEMORY]: $free_mem"
echo "Running: genomad end-to-end"

# Parse parameters:
#   -p  Input sequence file (e.g., .fna or .fna.gz)
#   -o  Output directory
#   -s  Number of splits (optional, default: 8)
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -p)
            input_file="$2"
            shift 2
            ;;
        -o)
            output_dir="$2"
            shift 2
            ;;
        -s)
            splits="$2"
            shift 2
            ;;
        *)
            echo "Unknown parameter: $1"
            exit 1
            ;;
    esac
done

# Check required parameters
if [ -z "$input_file" ]; then
    echo "Error: Input sequence file (-p) not provided."
    exit 1
fi

if [ -z "$output_dir" ]; then
    echo "Error: Output directory (-o) not provided."
    exit 1
fi

# If splits is not specified, default to 8
if [ -z "$splits" ]; then
    splits=8
fi

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Activate the genomad conda environment (ensure the environment with genomad is correctly installed)
echo "[INFO] Activating genomad environment 'genomad_env'..."
eval "$(conda shell.bash hook)"
conda activate genomad_env
if [ $? -ne 0 ]; then
    echo "Error: Failed to activate conda environment 'genomad_env'."
    exit 1
fi

# Check and download the database (default location: 'genomad_db' folder in the current directory)
db_dir="genomad_db"
if [ ! -d "$db_dir" ]; then
    echo "[INFO] Downloading genomad database to '$db_dir'..."
    genomad download-database "$db_dir"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to download the genomad database."
        exit 1
    fi
fi

# Execute the genomad analysis
echo "[INFO] Running genomad end-to-end (splits=$splits)..."
command="genomad end-to-end --cleanup --splits $splits \"$input_file\" \"$output_dir\" \"$db_dir\""
echo "[INFO] Executing command: $command"
eval $command
if [ $? -ne 0 ]; then
    echo "Error: genomad analysis pipeline failed."
    exit 1
fi

echo "genomad analysis pipeline completed successfully."