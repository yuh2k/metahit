# Generate Index
bwa index output/assembly/final_assembly.fasta

./bin/metahit-modules/alignment.sh \
    -f output/assembly/final_assembly.fasta \
    -1 output/readqc/hic/final_pure_reads_1.fastq \
    -2 output/readqc/hic/final_pure_reads_1.fastq \
    -o output/alignment
