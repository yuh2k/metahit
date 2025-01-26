#!/usr/bin/env bash
echo "[FREE MEMORY]: $free_mem"
echo "[INFO] Running Bin Plot"

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --contact-map) contact_map=$2; shift ;;
        --BIN) bin_file=$2; shift ;;
        --OUTDIR) out_dir=$2; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Check for required arguments
if [ -z "$contact_map" ] || [ -z "$bin_file" ] || [ -z "$out_dir" ]; then
    echo "[ERROR] Missing required arguments."
    echo "Usage: bin_plot.sh --contact-map <CONTACT_MAP> --BIN <BIN_FILE> --OUTDIR <OUTPUT_DIR>"
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

# Run the Python bin plot module
python "${script_dir}/bin/metahit-modules/bin_plot.py" \
    --contact-map "$contact_map" \
    --BIN "$bin_file" \
    --OUTDIR "$out_dir"

if [ $? -ne 0 ]; then
    echo "[ERROR] Bin Plot step failed."
    exit 1
fi

echo "[INFO] Bin Plot completed successfully."
