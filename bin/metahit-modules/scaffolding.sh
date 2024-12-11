#!/usr/bin/env bash

# If no arguments, show usage
if [ "$#" -lt 1 ]; then
    echo " -p metahit path"
    echo "Usage: $0 [options]"
    echo ""
    exit 1
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -p)
            path=$2
            shift 2
            ;;
        *)
            # Stop parsing here
            break
            ;;
    esac
done
mkdir "${path}/output/scaffolding/yahs"
touch final.bin
SCAFFOLDING_SCRIPT="${path}/bin/metahit-scripts/scaffolding.py"
if [ ! -f "$SCAFFOLDING_SCRIPT" ]; then
    echo "Error: scaffolding.py not found at $SCAFFOLDING_SCRIPT"
    exit 1
fi

python "$SCAFFOLDING_SCRIPT" "$@"
if [ $? -ne 0 ]; then
    echo "Error: Scaffolding step failed."
    exit 1
fi

echo "Scaffolding step completed successfully."
