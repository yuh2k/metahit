#!/usr/bin/env bash

help_message () {
    echo ""
    echo "Usage: metahit assembly [options] -1 reads_1.fastq -2 reads_2.fastq -o output_dir"
    echo "Options:"
    echo "    -1 STR          Forward fastq reads"
    echo "    -2 STR          Reverse fastq reads"
    echo "    -o STR          Output directory"
    echo "    -m INT          Memory in GB (default=24)" ####Todo available 80% available memory
    echo "    -t INT          Number of threads (default=1)" ####Todo 20
    echo "    -l INT          Minimum length of assembled contigs (default=1000)"
    echo ""
    echo "    --megahit       Assemble with MEGAHIT (default)"
    echo "    --k-min INT     MEGAHIT minimum k-mer size (default=21)"
    echo "    --k-max INT     MEGAHIT maximum k-mer size (default=141)"
    echo "    --k-step INT    MEGAHIT k-mer step size (default=12)"
    echo ""
    echo "    --metaspades    Assemble with metaSPAdes instead of MEGAHIT"
    echo "    --k-list STR    metaSPAdes k-mer sizes (e.g., '21,33,55,61'; default='21,33,55,61')"
    echo ""
}

SOFT="bin/metahit-scripts"

# Default parameters
mem=24
threads=1
out="false"
reads_1="false"
reads_2="false"
min_len=1000

metaspades_assemble=false
megahit_assemble=true

# MEGAHIT default parameters
k_min=21
k_max=141
k_step=12

# MetaSPAdes default parameters
k_list="21,33,55,61"

# Load in parameters
OPTS=$(getopt -o ht:m:o:1:2:l: --long help,megahit,k-min:,k-max:,k-step:,metaspades,k-list: -- "$@")
if [ $? != 0 ]; then help_message; exit 1; fi

eval set -- "$OPTS"

while true; do
    case "$1" in
        -t) threads=$2; shift 2;;
        -m) mem=$2; shift 2;;
        -o) out=$2; shift 2;;
        -1) reads_1=$2; shift 2;;
        -2) reads_2=$2; shift 2;;
        -l) min_len=$2; shift 2;;
        -h | --help) help_message; exit 0; shift ;;
        --megahit) megahit_assemble=true; metaspades_assemble=false; shift ;;
        --metaspades) metaspades_assemble=true; megahit_assemble=false; shift ;;
        --k-min) k_min=$2; shift 2;;
        --k-max) k_max=$2; shift 2;;
        --k-step) k_step=$2; shift 2;;
        --k-list) k_list=$2; shift 2;;
        --) shift; break ;;
        *) echo "Unknown option: $1"; help_message; exit 1;;
    esac
done

# Validate required arguments
if [ "$out" = "false" ] || [ "$reads_1" = "false" ] || [ "$reads_2" = "false" ]; then 
    help_message; exit 1
fi

# Dependency checks
echo "[INFO] Checking dependencies..."
if [ "$metaspades_assemble" = true ] && ! command -v metaspades.py &> /dev/null; then
    echo "[ERROR] metaSPAdes is not installed or not in the PATH."; exit 1
fi

if [ "$megahit_assemble" = true ] && ! command -v megahit &> /dev/null; then
    echo "[ERROR] MEGAHIT is not installed or not in the PATH."; exit 1
fi

if ! command -v quast.py &> /dev/null; then
    echo "[ERROR] QUAST is not installed or not in the PATH."; exit 1
fi

if [ ! -f "${SOFT}/rm_short_contigs.py" ]; then
    echo "[ERROR] Script 'rm_short_contigs.py' not found in ${SOFT}."; exit 1
fi

mkdir -p "$out" || { echo "Error: Cannot create output directory."; exit 1; }

# Validate and set up assembler-specific parameters
if [ "$metaspades_assemble" = true ]; then
    # Validate k-mer sizes for metaSPAdes
    IFS=',' read -ra K_ARRAY <<< "$k_list"
    for k in "${K_ARRAY[@]}"; do
        if ! [[ $k =~ ^[0-9]+$ ]]; then
            echo "Error: k-mer size '$k' is not an integer."
            exit 1
        fi
        if (( k % 2 == 0 )); then
            echo "Error: k-mer size '$k' is not odd."
            exit 1
        fi
        if (( k >= 128 )); then
            echo "Error: k-mer size '$k' must be less than 128."
            exit 1
        fi
    done
fi

if [ "$megahit_assemble" = true ]; then
    # Validate k-mer sizes for MEGAHIT
    if ! [[ $k_min =~ ^[0-9]+$ ]] || ! [[ $k_max =~ ^[0-9]+$ ]] || ! [[ $k_step =~ ^[0-9]+$ ]]; then
        echo "Error: k-min, k-max, and k-step must be integers."
        exit 1
    fi
    if (( k_min % 2 == 0 )); then
        echo "Error: MEGAHIT k-min '$k_min' is not odd."
        exit 1
    fi
    if (( k_max % 2 == 0 )); then
        echo "Error: MEGAHIT k-max '$k_max' is not odd."
        exit 1
    fi
    if (( k_min < 15 )); then
        echo "Error: MEGAHIT k-min '$k_min' must be at least 15."
        exit 1
    fi
    if (( k_max >= 255 )); then
        echo "Error: MEGAHIT k-max '$k_max' must be less than 255."
        exit 1
    fi
    if (( k_min > k_max )); then
        echo "Error: MEGAHIT k-min '$k_min' cannot be greater than k-max '$k_max'."
        exit 1
    fi
fi

# Run the selected assembler
if [ "$metaspades_assemble" = true ]; then
    echo "ASSEMBLING WITH metaSPAdes"
    mkdir -p "${out}/metaspades.tmp"
    echo "Using k-mer sizes: $k_list"
    metaspades.py --tmp-dir "${out}/metaspades.tmp" -k "$k_list" -t "$threads" -m "$mem" -o "${out}/metaspades" -1 "$reads_1" -2 "$reads_2"
    if [ -d "${out}/metaspades.tmp" ]; then rm -r "${out}/metaspades.tmp"; fi
    if [ ! -f "${out}/metaspades/scaffolds.fasta" ]; then echo "Error: metaSPAdes assembly failed."; exit 1; fi
fi

if [ "$megahit_assemble" = true ]; then
    echo "ASSEMBLING WITH MEGAHIT"
    echo "Using k-mer sizes from $k_min to $k_max with step size $k_step"
    megahit -1 "$reads_1" -2 "$reads_2" -o "${out}/megahit" --min-contig-len "$min_len" --k-min "$k_min" --k-max "$k_max" --k-step "$k_step" --merge-level 20,0.95 -t "$threads" -m "${mem}000000000"
    if [ ! -f "${out}/megahit/final.contigs.fa" ]; then echo "Error: MEGAHIT assembly failed."; exit 1; fi ## Todo: merge level 
fi

# Process the assembly output
if [ "$metaspades_assemble" = true ]; then
    ${SOFT}/rm_short_contigs.py "$min_len" "${out}/metaspades/scaffolds.fasta" > "${out}/final_assembly.fasta"
fi

if [[ ! -s "${out}/final_assembly.fasta" ]]; then echo "Error: Final assembly failed."; exit 1; fi

echo "RUNNING ASSEMBLY QC WITH QUAST"
quast.py -t "$threads" -o "${out}/QUAST_out" -m 500 "${out}/final_assembly.fasta" ## Todo -m 500
cp "${out}/QUAST ## Todo Check this command
