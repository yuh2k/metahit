#!/usr/bin/env bash

if [ "$#" -lt 1 ]; then
    echo " -p metahit path"
    echo "Usage: $0 <command> [options]"
    echo "Available commands:"
    echo "  run            Perform reassembly"
    echo ""
    echo "For help on each command, use:"
    echo "  $0 <command> --help"
    exit 1
fi

COMMAND=$1
shift

# Parse additional options
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

REASSEMBLY_SCRIPT="./${path}/bin/metahit-scripts/reassembly.py"
if [ ! -f "$REASSEMBLY_SCRIPT" ]; then
    echo "Error: reassembly.py not found at $REASSEMBLY_SCRIPT"
    exit 1
fi

python "$REASSEMBLY_SCRIPT" "$COMMAND" "$@"
if [ $? -ne 0 ]; then
    echo "Error: Reassembly step '$COMMAND' failed."
    exit 1
fi

echo "Reassembly step '$COMMAND' completed successfully."
