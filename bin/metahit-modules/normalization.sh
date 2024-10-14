#!/usr/bin/env bash

#############################################################################################################
# Normalization script using the normalization.py Python module.
#
# This script wraps the Python Normalization class, allowing normalization steps to be executed
# via command-line commands within the MetaHit pipeline.
#############################################################################################################

# Display usage if not enough arguments
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <command> [options]"
    echo "Available commands:"
    echo "  raw            Perform raw normalization"
    echo "  normcc_contact Perform normcc_contact normalization"
    echo "  hiczin         Perform HiCzin normalization"
    echo "  bin3c          Perform bin3C normalization"
    echo "  metator        Perform MetaTOR normalization"
    echo "  denoise        Perform denoising on a contact matrix"
    echo ""
    echo "For help on each command, use:"
    echo "  $0 <command> --help"
    exit 1
fi

COMMAND=$1
shift 1

# Path to the normalization.py script
NORMALIZATION_SCRIPT="./bin/metahit-scripts/normalization.py"

# Check if the normalization script exists
if [ ! -f "$NORMALIZATION_SCRIPT" ]; then
    echo "Error: Normalization script not found at $NORMALIZATION_SCRIPT"
    exit 1
fi

# Activate Conda environment (adjust the path if necessary)
source ~/anaconda3/etc/profile.d/conda.sh
conda activate metahit_env

# Execute the corresponding Python command
python "$NORMALIZATION_SCRIPT" "$COMMAND" "$@"

# Check exit status
if [ $? -ne 0 ]; then
    echo "Error: Normalization step '$COMMAND' failed."
    exit 1
fi

echo "Normalization step '$COMMAND' completed successfully."
