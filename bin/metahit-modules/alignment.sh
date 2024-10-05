"""
line 27: bwa: command not found
Samtools: -F 0x904 become option
make a new folder named alignment
"""
#!/usr/bin/env bash

#############################################################################################################
# Alignment script using BWA MEM and SamTools.
#
# Aligns reads to the final assembly and processes the output using SamTools.
#############################################################################################################

# Display usage if not enough arguments
if [ "$#" -ne 8 ]; then
    echo "Usage: $0 -f final_assembly.fasta -1 HIC_1_dedup.fastq.gz -2 HIC_2_dedup.fastq.gz -o output_dir"
    exit 1
fi

# Read in the arguments
while getopts f:1:2:o: flag
do
    case "${flag}" in
        f) fasta=${OPTARG};;
        1) reads_1=${OPTARG};;
        2) reads_2=${OPTARG};;
        o) output_dir=${OPTARG};;
    esac
done

mkdir -p "$output_dir"
bwa mem -5SP "$fasta" "$reads_1" "$reads_2" > "$output_dir/MAP.sam"
samtools view -F 0x904 -bS "$output_dir/MAP.sam" > "$output_dir/MAP_UNSORTED.bam"
samtools sort -n "$output_dir/MAP_UNSORTED.bam" -o "$output_dir/MAP_SORTED.bam"

echo "Alignment completed successfully!"
