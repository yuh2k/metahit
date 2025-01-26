#!/usr/bin/env bash
echo "[FREE MEMORY]: $free_mem"
if [ "$#" -lt 1 ]; then
    echo " -p metahit path"
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

COMMAND=$1
shift

# Parse additional options for this script
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -p) 
            path=$2
            shift 2
            ;;
        *)
            # Stop parsing here and pass the rest to python
            break
            ;;
    esac
done

# Path to the normalization.py script
NORMALIZATION_SCRIPT="${path}/bin/metahit-scripts/normalization.py"

if [ ! -f "$NORMALIZATION_SCRIPT" ]; then
    echo "Error: normalization.py not found at $NORMALIZATION_SCRIPT"
    exit 1
fi

# Pass all remaining arguments directly to python
python "$NORMALIZATION_SCRIPT" "$COMMAND" "$@"

if [ $? -ne 0 ]; then
    echo "Error: Normalization step '$COMMAND' failed."
    exit 1
fi

echo "Normalization step '$COMMAND' completed successfully."
