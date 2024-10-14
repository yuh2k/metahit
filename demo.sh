#!/usr/bin/env bash

#############################################################################################################
# MetaHit Workflow Demo Run Script
#############################################################################################################
# This script demonstrates a complete run of the MetaHit pipeline, including Read QC,
# Assembly, Alignment, Coverage Estimation, and Normalization. It accepts input parameters for files
# and directories, allowing flexible usage.
#############################################################################################################

set -e  # Exit immediately if a command exits with a non-zero status
set -o pipefail  # Ensure that pipeline errors are not masked

# Define input parameters
SG_R1="./test_data/sim_sg_R1.fastq"
SG_R2="./test_data/sim_sg_R2.fastq"
HIC_R1="./test_data/sim_hic_r1.fq"
HIC_R2="./test_data/sim_hic_r2.fq"
OUTPUT_DIR="output"

# Activate Conda environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate metahit_env

# ====================================================================================
# Read QC
# ====================================================================================

# Metagenomic shotgun sequencing
./bin/metahit-modules/read_qc.sh -1 "$SG_R1" -2 "$SG_R2" -o "${OUTPUT_DIR}/readqc/sg" -t 4

# Hi-C sequencing
./bin/metahit-modules/read_qc.sh -1 "$HIC_R1" -2 "$HIC_R2" -o "${OUTPUT_DIR}/readqc/hic" -t 4

# ====================================================================================
# Assembly
# ====================================================================================

./bin/metahit-modules/assembly.sh \
    -1 "${OUTPUT_DIR}/readqc/sg/final_pure_reads_1.fastq" \
    -2 "${OUTPUT_DIR}/readqc/sg/final_pure_reads_2.fastq" \
    -o "${OUTPUT_DIR}/assembly" \
    -m 24 -t 4 --megahit

# ====================================================================================
# Generate Index
# ====================================================================================

bwa index "${OUTPUT_DIR}/assembly/final_assembly.fasta"

# ====================================================================================
# Alignment
# ====================================================================================

./bin/metahit-modules/alignment.sh \
    -f "${OUTPUT_DIR}/assembly/final_assembly.fasta" \
    -1 "${OUTPUT_DIR}/readqc/hic/final_pure_reads_1.fastq" \
    -2 "${OUTPUT_DIR}/readqc/hic/final_pure_reads_2.fastq" \
    -o "${OUTPUT_DIR}/alignment"

# ====================================================================================
# Coverage Estimation
# ====================================================================================

./bin/metahit-modules/coverage_estimation.sh \
    -1 "$SG_R1" \
    -2 "$SG_R2" \
    -r "${OUTPUT_DIR}/assembly/final_assembly.fasta" \
    -o "${OUTPUT_DIR}/estimation"

# ====================================================================================
# Normalization Steps
# ====================================================================================

# Define normalization parameters
BAM_FILE="./output/alignment/sorted_map.bam"
FASTA_FILE="./output/assembly/final_assembly.fasta"
ENZYMES=("EcoRI" "HindIII")
NORMALIZATION_OUTPUT="./output/normalization"

# 1. Raw Normalization
./bin/metahit-modules/normalization.sh raw \
    -b "$BAM_FILE" \
    -f "$FASTA_FILE" \
    -e "${ENZYMES[@]}" \
    -o "${NORMALIZATION_OUTPUT}/raw" \
    --min_mapq 20 \
    --min_len 500 \
    --min_match 15 \
    --min_signal 1

# 2. normcc_contact Normalization
./bin/metahit-modules/normalization.sh normcc_contact \
    -c "${NORMALIZATION_OUTPUT}/raw/contig_info.csv" \
    -m "${NORMALIZATION_OUTPUT}/raw/contact_matrix.npz" \
    -o "${NORMALIZATION_OUTPUT}/normcc"

# 3. HiCzin Normalization
./bin/metahit-modules/normalization.sh hiczin \
    -c "${NORMALIZATION_OUTPUT}/raw/contig_info.csv" \
    -m "${NORMALIZATION_OUTPUT}/raw/contact_matrix.npz" \
    -o "${NORMALIZATION_OUTPUT}/hiczin" \
    --min_signal 1

# 4. bin3C Normalization
./bin/metahit-modules/normalization.sh bin3c \
    -c "${NORMALIZATION_OUTPUT}/raw/contig_info.csv" \
    -m "${NORMALIZATION_OUTPUT}/raw/contact_matrix.npz" \
    -o "${NORMALIZATION_OUTPUT}/bin3c" \
    --max_iter 1000 \
    --tol 1e-6

# 5. MetaTOR Normalization
./bin/metahit-modules/normalization.sh metator \
    -c "${NORMALIZATION_OUTPUT}/raw/contig_info.csv" \
    -m "${NORMALIZATION_OUTPUT}/raw/contact_matrix.npz" \
    -o "${NORMALIZATION_OUTPUT}/metator"

# ====================================================================================
# Denoising Steps
# ====================================================================================

# Define normalization methods and their contact matrices
declare -A normalization_methods=(
    ["raw"]="${NORMALIZATION_OUTPUT}/raw/contact_matrix.npz"
    ["normcc"]="${NORMALIZATION_OUTPUT}/normcc/normalized_contact_matrix.npz"
    ["hiczin"]="${NORMALIZATION_OUTPUT}/hiczin/hiczin_contact_matrix.npz"
    ["bin3c"]="${NORMALIZATION_OUTPUT}/bin3c/bin3c_contact_matrix.npz"
    ["metator"]="${NORMALIZATION_OUTPUT}/metator/metator_contact_matrix.npz"
)

# Apply denoising to each normalization method
for method in "${!normalization_methods[@]}"; do
    cm_file=${normalization_methods[$method]}
    denoise_path="${NORMALIZATION_OUTPUT}/${method}/denoise"

    if [[ -f "$cm_file" ]]; then
        echo "Applying denoising to ${method} normalization."
        ./bin/metahit-modules/normalization.sh denoise \
            -m "$cm_file" \
            -o "$denoise_path" \
            -p 0.1
        echo "Denoising for ${method} normalization completed successfully."
    else
        echo "Contact matrix file for ${method} normalization does not exist. Skipping denoising."
    fi
done

echo "MetaHit pipeline demo run completed successfully!"
