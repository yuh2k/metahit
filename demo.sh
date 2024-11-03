#!/bin/bash

# Activate Conda environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate metahit_env

# Define input parameters
SG_R1="./test_data/sim_sg_R1.fastq"
SG_R2="./test_data/sim_sg_R2.fastq"
HIC_R1="./test_data/sim_hic_r1.fq"
HIC_R2="./test_data/sim_hic_r2.fq"
OUTPUT_DIR="output"

# Run each command with the new Python CLI

# Read QC
./metahit.py readqc -1 "$SG_R1" -2 "$SG_R2" -o "${OUTPUT_DIR}/readqc/sg" -t 4
./metahit.py readqc -1 "$HIC_R1" -2 "$HIC_R2" -o "${OUTPUT_DIR}/readqc/hic" -t 4

# Assembly
./metahit.py assembly -1 "${OUTPUT_DIR}/readqc/sg/final_reads_1.fastq" -2 "${OUTPUT_DIR}/readqc/sg/final_reads_2.fastq" -o "${OUTPUT_DIR}/assembly" -m 24 -t 4 --metaspades --k-list 21,33,55,77 -l 1000

# Alignment
./metahit.py alignment -r "${OUTPUT_DIR}/assembly/final_assembly.fasta" -1 "${OUTPUT_DIR}/readqc/hic/final_reads_1.fastq" -2 "${OUTPUT_DIR}/readqc/hic/final_reads_2.fastq" -o "${OUTPUT_DIR}/alignment" --threads 4 --samtools-filter '-F 0x904'

# Coverage Estimation
./metahit.py coverage_estimation -1 "$SG_R1" -2 "$SG_R2" -r "${OUTPUT_DIR}/assembly/final_assembly.fasta" -o "${OUTPUT_DIR}/estimation"

# Raw Contact Generation
./metahit.py raw_contact --bam "${OUTPUT_DIR}/alignment/sorted_map.bam" --fasta "${OUTPUT_DIR}/assembly/final_assembly.fasta" --out "${OUTPUT_DIR}/normalization/raw" --coverage "${OUTPUT_DIR}/estimation/coverage.txt"

# Normalization
./metahit.py normalization raw --contig_file "${OUTPUT_DIR}/normalization/raw/contig_info.csv" --contact_matrix_file "${OUTPUT_DIR}/normalization/raw/contact_matrix_user.npz" --output "${OUTPUT_DIR}/normalization/raw_normalized" --min_len 500 --min_signal 1 --thres 5
./metahit.py normalization normcc --contig_file "${OUTPUT_DIR}/normalization/raw/contig_info.csv" --contact_matrix_file "${OUTPUT_DIR}/normalization/raw/contact_matrix_user.npz" --output "${OUTPUT_DIR}/normalization/normcc" --thres 5
./metahit.py normalization hiczin --contig_file "${OUTPUT_DIR}/normalization/raw/contig_info.csv" --contact_matrix_file "${OUTPUT_DIR}/normalization/raw/contact_matrix_user.npz" --output "${OUTPUT_DIR}/normalization/hiczin" --thres 5
./metahit.py normalization bin3c --contig_file "${OUTPUT_DIR}/normalization/raw/contig_info.csv" --contact_matrix_file "${OUTPUT_DIR}/normalization/raw/contact_matrix_user.npz" --output "${OUTPUT_DIR}/normalization/bin3c" --max_iter 1000 --tol 1e-6 --thres 5
./metahit.py normalization metator --contig_file "${OUTPUT_DIR}/normalization/raw/contig_info.csv" --contact_matrix_file "${OUTPUT_DIR}/normalization/raw/contact_matrix_user.npz" --output "${OUTPUT_DIR}/normalization/metator" --thres 5

