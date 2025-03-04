#!/usr/bin/env bash

###########################################################################################################
# This script is a comprehensive solution to QC HiSeq reads in preparation for assembly and other
# operations. The main tasks are read trimming based on quality scores and producing FASTQC reports.
###########################################################################################################

# Function to display the help message
help_message () {
    echo ""
    echo "Usage: read_qc.sh [options] -1 reads_1.fastq -2 reads_2.fastq -o output_dir"
    echo "Options:"
    echo "    -p metahit path"
    echo "    -1 STR          forward fastq reads (.fastq or .fastq.gz)"
    echo "    -2 STR          reverse fastq reads (.fastq or .fastq.gz)"
    echo "    -o STR          output directory"
    echo "    -t INT          number of threads (default=1)"
    echo "    --minlen INT    minimum length of reads after trimming (default=50)"
    echo "    --trimq INT     quality threshold for trimming (default=10)"
    echo "    --ftl INT       trim bases from the left (default=10)"
    echo "    --xmx STR       memory for Java (default=60% of available memory)"
    echo "    --ftm INT       optional ftm value (default=5)"
    echo "    --dedup         perform deduplication with clumpify (default=false)"
    echo "    --skip-pre-qc-report   skip FastQC for input sequences"
    echo "    --skip-post-qc-report  skip FastQC for final sequences"
    echo ""
}

# Defaults
threads=1
minlen=50
trimq=10
ftl=10
dedup=false
pre_qc_report=true
post_qc_report=true
out=""
reads_1=""
reads_2=""
ftm=5

# Calculate default xmx as 60% of available memory
available_mem_kb=$(grep MemTotal /proc/meminfo | awk '{print $2}')
default_xmx=$((available_mem_kb * 60 / 100 / 1024))g
xmx=$default_xmx

# Getopts
OPTS=$(getopt -o hp:t:o:1:2: --long help,skip-pre-qc-report,skip-post-qc-report,dedup,minlen:,trimq:,ftl:,xmx:,ftm: -- "$@")
if [ $? -ne 0 ]; then help_message; exit 1; fi
eval set -- "$OPTS"

# Parse options
while true; do
    case "$1" in
        -p) path=$2; shift 2;;
        -t) threads=$2; shift 2;;
        -o) out=$2; shift 2;;
        -1) reads_1=$2; shift 2;;
        -2) reads_2=$2; shift 2;;
        --minlen) minlen=$2; shift 2;;
        --trimq) trimq=$2; shift 2;;
        --ftl) ftl=$2; shift 2;;
        --xmx) xmx=$2; shift 2;;
        --ftm) ftm=$2; shift 2;;
        --dedup) dedup=true; shift 1;;
        --skip-pre-qc-report) pre_qc_report=false; shift 1;;
        --skip-post-qc-report) post_qc_report=false; shift 1;;
        -h|--help) help_message; exit 0;;
        --) shift; break;;
        *) echo "Invalid option: $1"; help_message; exit 1;;
    esac
done

# Check if all required parameters are entered
if [ -z "$out" ] || [ -z "$reads_1" ] || [ -z "$reads_2" ]; then
    help_message; exit 1
fi

# Check if output directory exists, if not, create it
mkdir -p "$out"

# Check if input read files exist in either .fastq or .fastq.gz formats
for read_file in "$reads_1" "$reads_2"; do
    if [[ ! -f "$read_file" ]] && [[ ! -f "${read_file}.gz" ]]; then
        echo "Error: $read_file or ${read_file}.gz not found"
        exit 1
    fi
done

# FastQC report for input reads if requested
if [ "$pre_qc_report" = true ]; then
    mkdir -p "${out}/pre-QC_report"
    fastqc -q -t "$threads" -o "${out}/pre-QC_report" "$reads_1" "$reads_2"
    rm "${out}/pre-QC_report"/*zip
fi

ref="${path}/external/bbmap/resources/adapters.fa"

# Run BBDuk Steps with .gz output
echo "Running Adapter trimming with BBDuk"
${path}/external/bbmap/bbduk.sh -Xmx$xmx in1="$reads_1" in2="$reads_2" out1="${out}/step1_adptrim_1.fastq.gz" out2="${out}/step1_adptrim_2.fastq.gz" \
 ref="$ref" ktrim=r k=23 mink=11 hdist=1 minlen="$minlen" threads="$threads" gzip=true

echo "Running Quality trimming with BBDuk"
${path}/external/bbmap/bbduk.sh -Xmx$xmx in1="${out}/step1_adptrim_1.fastq.gz" in2="${out}/step1_adptrim_2.fastq.gz" \
    out1="${out}/step2_qualtrim_1.fastq.gz" out2="${out}/step2_qualtrim_2.fastq.gz" qtrim=r trimq="$trimq" ftm="$ftm" \
    minlen="$minlen" threads="$threads" gzip=true

echo "Trimming left bases with BBDuk"
${path}/external/bbmap/bbduk.sh -Xmx$xmx in1="${out}/step2_qualtrim_1.fastq.gz" in2="${out}/step2_qualtrim_2.fastq.gz" \
    out1="${out}/step3_lefttrim_1.fastq.gz" out2="${out}/step3_lefttrim_2.fastq.gz" ftl="$ftl" threads="$threads" gzip=true

# Deduplication step
if [ "$dedup" = true ]; then
    echo "Running deduplication with Clumpify"
    clumpify.sh in1="${out}/step3_lefttrim_1.fastq.gz" in2="${out}/step3_lefttrim_2.fastq.gz" \
        out1="${out}/final_reads_1.fastq.gz" out2="${out}/final_reads_2.fastq.gz" dedupe threads="$threads" gzip=true
else
    mv "${out}/step3_lefttrim_1.fastq.gz" "${out}/final_reads_1.fastq.gz"
    mv "${out}/step3_lefttrim_2.fastq.gz" "${out}/final_reads_2.fastq.gz"
fi

# FastQC report for final reads if requested
if [ "$post_qc_report" = true ]; then
    mkdir -p "${out}/post-QC_report"
    fastqc -t "$threads" -o "${out}/post-QC_report" "${out}/final_reads_1.fastq.gz" "${out}/final_reads_2.fastq.gz"
    rm "${out}/post-QC_report"/*zip
fi

echo "READ QC PIPELINE COMPLETE!"
