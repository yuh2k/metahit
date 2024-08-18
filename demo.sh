# Read QC

#  metagenomic shotgun sequencing 
./bin/metahit-modules/read_qc.sh -1 test_data/sim_sg_R1.fastq -2 test_data/sim_sg_R2.fastq -o output/readqc/sg -t 4
#  Hi-C sequencing 
./bin/metahit-modules/read_qc.sh -1 test_data/sim_hic_r1.fq -2 test_data/sim_hic_r2.fq -o output/readqc/hic -t 4

# =======================================================================================================================

# Assembly
./bin/metahit-modules/assembly.sh -1 output/readqc/sg/final_pure_reads_1.fastq -2 output/readqc/sg/final_pure_reads_2.fastq -o output/assembly -m 24 -t 4 --megahit

# =======================================================================================================================

# Generate Index
bwa index output/assembly/final_assembly.fasta

./bin/metahit-modules/alignment.sh \
    -f output/assembly/final_assembly.fasta \
    -1 output/readqc/hic/final_pure_reads_1.fastq \
    -2 output/readqc/hic/final_pure_reads_1.fastq \
    -o output/alignment

# =======================================================================================================================

./bin/metahit-modules/coverage_estimation.sh -1 test_data/sim_sg_R1.fastq -2 test_data/sim_sg_R2.fastq -r output/assembly/final_assembly.fasta -o output/estimation
