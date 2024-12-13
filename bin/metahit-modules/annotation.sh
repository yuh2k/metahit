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

# Activate conda environment
eval "$(conda shell.bash hook)"
if ! conda activate gtdbtk-2.1.1; then
    echo "Error: Could not activate conda environment 'gtdbtk-2.1.1'."
    exit 1
fi

# Run GTDB-Tk annotation
echo "Running: gtdbtk classify_wf $@"
gtdbtk classify_wf "$@"
if [ $? -ne 0 ]; then
    echo "Error: Annotation (GTDB-Tk) step failed."
    exit 1
fi

echo "Annotation (GTDB-Tk) step completed successfully."
