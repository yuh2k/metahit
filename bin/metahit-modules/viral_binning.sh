#!/usr/bin/env bash

if [ "$#" -lt 1 ]; then
    echo " -p metahit path"
    echo "Usage: $0 <command> [options]"
    echo "Available commands:"
    echo "  pipeline       Run the viralcc pipeline"
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
            break
            ;;
    esac
done
VIRALCC_SCRIPT="${path}/bin/metahit-scripts/viralcc/viralcc.py"
if [ ! -f "$VIRALCC_SCRIPT" ]; then
    echo "Error: viralcc.py not found at $VIRALCC_SCRIPT"
    exit 1
fi

python "$VIRALCC_SCRIPT" "$COMMAND" "$@"
if [ $? -ne 0 ]; then
    echo "Error: ViralCC step '$COMMAND' failed."
    exit 1
fi

echo "ViralCC step '$COMMAND' completed successfully."
