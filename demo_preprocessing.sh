#!/usr/bin/env bash

# ====================================================================================
# MetaHit Workflow Demo Run Script
# ====================================================================================
# This script demonstrates a complete run of the MetaHit pipeline, including Read QC,
# Assembly, Alignment, and Coverage Estimation. It accepts input parameters for files
# and directories, allowing flexible usage.
# ====================================================================================

# Define input parameters
SG_R1="test_data/sim_sg_R1.fastq"
SG_R2="test_data/sim_sg_R2.fastq"
HIC_R1="test_data/sim_hic_r1.fq"
HIC_R2="test_data/sim_hic_r2.fq"
OUTPUT_DIR="output"

# ====================================================================================
# Read QC
# ====================================================================================

# Metagenomic shotgun sequencing
./bin/metahit-modules/read_qc.sh -1 $SG_R1 -2 $SG_R2 -o ${OUTPUT_DIR}/readqc/sg -t 4

# Hi-C sequencing
./bin/metahit-modules/read_qc.sh -1 $HIC_R1 -2 $HIC_R2 -o ${OUTPUT_DIR}/readqc/hic -t 4

# ====================================================================================
# Assembly
# ====================================================================================

./bin/metahit-modules/assembly.sh \
    -1 ${OUTPUT_DIR}/readqc/sg/final_pure_reads_1.fastq \
    -2 ${OUTPUT_DIR}/readqc/sg/final_pure_reads_2.fastq \
    -o ${OUTPUT_DIR}/assembly \
    -m 24 -t 4 --megahit

# ====================================================================================
# Generate Index
# ====================================================================================

bwa index ${OUTPUT_DIR}/assembly/final_assembly.fasta

# ====================================================================================
# Alignment
# ====================================================================================

./bin/metahit-modules/alignment.sh \
    -f ${OUTPUT_DIR}/assembly/final_assembly.fasta \
    -1 ${OUTPUT_DIR}/readqc/hic/final_pure_reads_1.fastq \
    -2 ${OUTPUT_DIR}/readqc/hic/final_pure_reads_2.fastq \
    -o ${OUTPUT_DIR}/alignment

# ====================================================================================
# Coverage Estimation
# ====================================================================================

./bin/metahit-modules/coverage_estimation.sh \
    -1 $SG_R1 \
    -2 $SG_R2 \
    -r ${OUTPUT_DIR}/assembly/final_assembly.fasta \
    -o ${OUTPUT_DIR}/estimation

echo "MetaHit pipeline demo run completed successfully!"
