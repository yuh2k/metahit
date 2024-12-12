#!/usr/bin/env bash

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
echo "Running GTDB"
mkdir 123123
mkdir ./{$path}/output/annotation
eval "$(conda shell.bash hook)"
conda activate gtdbtk-2.1.0

gtdbtk classify_wf "$@"
if [ $? -ne 0 ]; then
    echo "Error: Annotation (GTDB-Tk) step failed."
    exit 1
fi

echo "Annotation (GTDB-Tk) step completed successfully."
