#!/usr/bin/env bash

# Activate Conda environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate metahit_env

# Define input parameters
SG_R1="./test_data/sim_sg_R1.fastq"
SG_R2="./test_data/sim_sg_R2.fastq"
HIC_R1="./test_data/sim_hic_r1.fq"
HIC_R2="./test_data/sim_hic_r2.fq"
OUTPUT_DIR="output"
REFINEMENT_METHOD="metacc"

# Step 1: Read QC
echo "[INFO] Running Read QC..."
./metahit.py readqc -1 "$SG_R1" -2 "$SG_R2" -o "${OUTPUT_DIR}/readqc/sg" -t 4 --xmx 4g
./metahit.py readqc -1 "$HIC_R1" -2 "$HIC_R2" -o "${OUTPUT_DIR}/readqc/hic" -t 4 --xmx 4g

# Step 2: Assembly
echo "[INFO] Running Assembly..."
./metahit.py assembly -1 "${OUTPUT_DIR}/readqc/sg/final_reads_1.fastq.gz" -2 "${OUTPUT_DIR}/readqc/sg/final_reads_2.fastq.gz" \
  -o "${OUTPUT_DIR}/assembly" -m 24 -t 4 --metaspades --k-list 21,33,55,77 -l 1000

# Step 3: Alignment
echo "[INFO] Running Alignment..."
./metahit.py alignment -r "${OUTPUT_DIR}/assembly/final_assembly.fasta" \
  -1 "${OUTPUT_DIR}/readqc/hic/final_reads_1.fastq.gz" -2 "${OUTPUT_DIR}/readqc/hic/final_reads_2.fastq.gz" \
  -o "${OUTPUT_DIR}/alignment" --threads 4 --samtools-filter '-F 0x904'

# Step 4: Coverage Estimation
echo "[INFO] Running Coverage Estimation..."
./metahit.py coverage_estimation -1 "$SG_R1" -2 "$SG_R2" \
  -r "${OUTPUT_DIR}/assembly/final_assembly.fasta" -o "${OUTPUT_DIR}/estimation"

# Step 5: Raw Contact Generation
echo "[INFO] Generating Raw Contacts..."
./metahit.py raw_contact \
    --bam "${OUTPUT_DIR}/alignment/sorted_map.bam" \
    --fasta "${OUTPUT_DIR}/assembly/final_assembly.fasta" \
    --out "${OUTPUT_DIR}/normalization/raw" \
    --enzyme "HindIII"

# Step 6: Normalization with optional refinement
echo "[INFO] Running Normalization with refinement..."
for method in raw normcc hiczin bin3c metator fastnorm; do
  output_dir="${OUTPUT_DIR}/normalization/${method}"

  # Ensure output directory exists
  mkdir -p "$output_dir"

  case $method in
    raw)
      ./metahit.py normalization raw --contig_file "${OUTPUT_DIR}/normalization/raw/contig_info.csv" \
        --contact_matrix_file "${OUTPUT_DIR}/normalization/raw/raw_contact_matrix_metacc.npz" \
        --output "${OUTPUT_DIR}/normalization/raw" --min_len 500 --min_signal 1 --thres 1
      ;;
    normcc)
      echo "[INFO] running normcc normalization"
      ./metahit.py normalization $method --contig_file "${OUTPUT_DIR}/normalization/raw/contig_info.csv" \
        --contact_matrix_file "${OUTPUT_DIR}/normalization/raw/raw_contact_matrix_metacc.npz" \
        --output "${OUTPUT_DIR}/normalization/${method}" --thres 1
      ;;
    hiczin)
      ./metahit.py normalization $method --contig_file "${OUTPUT_DIR}/normalization/raw/contig_info.csv" \
        --contact_matrix_file "${OUTPUT_DIR}/normalization/raw/raw_contact_matrix_metacc.npz" \
        --output "${OUTPUT_DIR}/normalization/${method}" --thres 1 --epsilon 1
      ;;
    metator)
      ./metahit.py normalization $method --contig_file "${OUTPUT_DIR}/normalization/raw/contig_info.csv" \
        --contact_matrix_file "${OUTPUT_DIR}/normalization/raw/raw_contact_matrix_metacc.npz" \
        --output "${OUTPUT_DIR}/normalization/${method}" --thres 1 --epsilon 1
      ;;
    bin3c)
      ./metahit.py normalization bin3c --contig_file "${OUTPUT_DIR}/normalization/raw/contig_info.csv" \
        --contact_matrix_file "${OUTPUT_DIR}/normalization/raw/raw_contact_matrix_metacc.npz" \
        --output "${OUTPUT_DIR}/normalization/bin3c" --max_iter 1000 --tol 1e-6 --thres 1
      ;;
    fastnorm)
      ./metahit.py normalization fastnorm --contig_file "${OUTPUT_DIR}/normalization/raw/contig_info.csv" \
        --contact_matrix_file "${OUTPUT_DIR}/normalization/raw/raw_contact_matrix_metacc.npz" \
        --output "${OUTPUT_DIR}/normalization/fastnorm" --epsilon 1 --thres 1
      ;;
  esac

  if [ $? -ne 0 ]; then
    echo "[ERROR] Normalization failed for method $method."
    exit 1
  fi
done

# Step 7: Binning
# echo "[INFO] Running Binning Process..."
# ./metahit.py binning --bam "${OUTPUT_DIR}/alignment/sorted_map.bam" \
#   --fasta "${OUTPUT_DIR}/assembly/final_assembly.fasta" \
#   --coverage "${OUTPUT_DIR}/estimation/coverage.txt" \
#   --output "${OUTPUT_DIR}/bins" \
#   --enzymes HindIII \
#   --num_gene 100 \
#   --seed 42

# echo "[INFO] All steps completed successfully."
