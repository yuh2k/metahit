#!/usr/bin/env bash

#############################################################################################################
# MetaHit Workflow Demo Run Script
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

# Adjusted assembly step to include new options for MEGAHIT
./bin/metahit-modules/assembly.sh \
    -1 "${OUTPUT_DIR}/readqc/sg/final_reads_1.fastq" \
    -2 "${OUTPUT_DIR}/readqc/sg/final_reads_2.fastq" \
    -o "${OUTPUT_DIR}/assembly" \
    -m 24 -t 4 \
    --megahit \
    --k-min 21 --k-max 141 --k-step 12 \
    -l 1000

# If you prefer to use metaSPAdes instead of MEGAHIT, uncomment the following lines and comment out the MEGAHIT section above
# ./bin/metahit-modules/assembly.sh \
#     -1 "${OUTPUT_DIR}/readqc/sg/final_reads_1.fastq" \
#     -2 "${OUTPUT_DIR}/readqc/sg/final_reads_2.fastq" \
#     -o "${OUTPUT_DIR}/assembly" \
#     -m 24 -t 4 \
#     --metaspades \
#     --k-list 21,33,55,77 \
#     -l 1000

# ====================================================================================
# Alignment
# ====================================================================================

# Adjusted alignment step to include the samtools filter option if needed
./bin/metahit-modules/alignment.sh \
    -r "${OUTPUT_DIR}/assembly/final_assembly.fasta" \
    -1 "${OUTPUT_DIR}/readqc/hic/final_reads_1.fastq" \
    -2 "${OUTPUT_DIR}/readqc/hic/final_reads_2.fastq" \
    -o "${OUTPUT_DIR}/alignment" \
    --threads 4 \
    --samtools-filter '-F 0x904'

# ====================================================================================
# Coverage Estimation
# ====================================================================================

./bin/metahit-modules/coverage_estimation.sh \
    -1 "$SG_R1" \
    -2 "$SG_R2" \
    -r "${OUTPUT_DIR}/assembly/final_assembly.fasta" \
    -o "${OUTPUT_DIR}/estimation"

# ====================================================================================
# Raw Contact Generation
# ====================================================================================

# Define normalization parameters
BAM_FILE="${OUTPUT_DIR}/alignment/sorted_map.bam"
FASTA_FILE="${OUTPUT_DIR}/assembly/final_assembly.fasta"
NORMALIZATION_OUTPUT="${OUTPUT_DIR}/normalization"
COVERAGE_FILE="${OUTPUT_DIR}/estimation/coverage.txt"

# Run raw_contact.sh to generate contact matrices and contig info
./bin/metahit-modules/raw_contact.sh \
    --bam "$BAM_FILE" \
    --fasta "$FASTA_FILE" \
    --out "${NORMALIZATION_OUTPUT}/raw" \
    --coverage "$COVERAGE_FILE"

# ====================================================================================
# Normalization Steps
# ====================================================================================

# Set common parameters
CONTIG_FILE="${NORMALIZATION_OUTPUT}/raw/contig_info.csv"
CONTACT_MATRIX_FILE="${NORMALIZATION_OUTPUT}/raw/contact_matrix_user.npz"

# Raw Normalization
./bin/metahit-modules/normalization.sh raw \
    --contig_file "$CONTIG_FILE" \
    --contact_matrix_file "$CONTACT_MATRIX_FILE" \
    --output_path "${NORMALIZATION_OUTPUT}/raw_normalized" \
    --min_len 500 \
    --min_signal 1 \
    --thres 5

# normCC Normalization
./bin/metahit-modules/normalization.sh normcc \
    --contig_file "$CONTIG_FILE" \
    --contact_matrix_file "$CONTACT_MATRIX_FILE" \
    --output_path "${NORMALIZATION_OUTPUT}/normcc" \
    --thres 5

# HiCzin Normalization
./bin/metahit-modules/normalization.sh hiczin \
    --contig_file "$CONTIG_FILE" \
    --contact_matrix_file "$CONTACT_MATRIX_FILE" \
    --output_path "${NORMALIZATION_OUTPUT}/hiczin" \
    --thres 5

# bin3C Normalization
./bin/metahit-modules/normalization.sh bin3c \
    --contig_file "$CONTIG_FILE" \
    --contact_matrix_file "$CONTACT_MATRIX_FILE" \
    --output_path "${NORMALIZATION_OUTPUT}/bin3c" \
    --max_iter 1000 \
    --tol 1e-6 \
    --thres 5

# MetaTOR Normalization
./bin/metahit-modules/normalization.sh metator \
    --contig_file "$CONTIG_FILE" \
    --contact_matrix_file "$CONTACT_MATRIX_FILE" \
    --output_path "${NORMALIZATION_OUTPUT}/metator" \
    --thres 5

echo "MetaHit pipeline demo run completed successfully!"
