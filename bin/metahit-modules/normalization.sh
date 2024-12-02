#!/usr/bin/env bash

# Display usage if not enough arguments
if [ "$#" -lt 1 ]; then
    echo "    -p metahit path"
    echo "Usage: $0 <command> [options]"
    echo "Available commands:"
    echo "  raw            Perform raw normalization"
    echo "  normcc         Perform normCC normalization"
    echo "  hiczin         Perform HiCzin normalization"
    echo "  bin3c          Perform bin3C normalization"
    echo "  metator        Perform MetaTOR normalization"
    echo ""
    echo "For help on each command, use:"
    echo "  $0 <command> --help"
    exit 1
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -p) path=$2; shift 2;;
    esac
    shift
done

COMMAND=$1
shift 1

# Path to the normalization.py script
NORMALIZATION_SCRIPT="${path}/bin/metahit-scripts/normalization.py"

# Execute the corresponding Python command
python "$NORMALIZATION_SCRIPT" "$COMMAND" "$@"

if [ $? -ne 0 ]; then
    echo "Error: Normalization step '$COMMAND' failed."
    exit 1
fi

echo "Normalization step '$COMMAND' completed successfully."
