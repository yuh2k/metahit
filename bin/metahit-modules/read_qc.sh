#!/usr/bin/env bash

###########################################################################################################
#                                                                                                         #
# This script is a comprehensive solution to QC new HiSeq reads in preparation for assembly and other     #
# operations. The main things this pipeline accomplishes are read trimming based on quality scores.       #
# The script also produces a FASTQC report before and after the procedures.                               #
#                                                                                                         #
###########################################################################################################

# Function to display the help message
help_message () {
    echo ""
    echo "Usage: read_qc.sh [options] -1 reads_1.fastq -2 reads_2.fastq -o output_dir"
    echo "Note: the read files have to be named in the name_1.fastq/name_2.fastq convention."
    echo "Options:"
    echo ""
    echo "    -1 STR          forward fastq reads"
    echo "    -2 STR          reverse fastq reads"
    echo "    -o STR          output directory"
    echo "    -t INT          number of threads (default=1)"
    echo "    --minlen INT    minimum length of reads after trimming (default=50)"
    echo "    --trimq INT     quality threshold for trimming (default=10)"
    echo "    --ftl INT       trim bases from the left (default=10)"
    echo "    --dedup         perform deduplication with clumpify (default=false)"
    echo "    --skip-pre-qc-report   don't make FastQC report of input sequences"
    echo "    --skip-post-qc-report  don't make FastQC report of final sequences"
    echo "";}

########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################

# Default parameters
threads=1
out=""
reads_1=""
reads_2=""
pre_qc_report=true
post_qc_report=true
dedup=false
minlen=50
trimq=10
ftl=10
ref="external/bbmap/resources/adapters.fa"

# Load in parameters using getopt
OPTS=$(getopt -o ht:o:1:2: --long help,skip-pre-qc-report,skip-post-qc-report,dedup,minlen:,trimq:,ftl: -- "$@")
# Make sure the parameters are entered correctly
if [ $? -ne 0 ]; then help_message; exit 1; fi

eval set -- "$OPTS"

# Loop through input parameters
while true; do
    case "$1" in
        -t) threads=$2; shift 2;;
        -o) out=$2; shift 2;;
        -1) reads_1=$2; shift 2;;
        -2) reads_2=$2; shift 2;;
        --minlen) minlen=$2; shift 2;;
        --trimq) trimq=$2; shift 2;;
        --ftl) ftl=$2; shift 2;;
        --dedup) dedup=true; shift 1;;
        --skip-pre-qc-report) pre_qc_report=false; shift 1;;
        --skip-post-qc-report) post_qc_report=false; shift 1;;
        -h|--help) help_message; exit 0;;
        --) shift; break;;
        *) echo "Invalid option: $1"; help_message; exit 1;;
    esac
done

########################################################################################################
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################

# Check if all required parameters are entered
if [ -z "$out" ] || [ -z "$reads_1" ] || [ -z "$reads_2" ]; then
    help_message; exit 1
fi

# Check if forward and reverse reads are the same file
if [ "$reads_1" = "$reads_2" ]; then
    echo "The forward and reverse reads are the same file. Exiting pipeline."
    exit 1
fi

# Check if output directory exists, if not, create it
if [ ! -d "$out" ]; then
    mkdir -p "$out"
else
    echo "Warning: $out already exists."
fi

########################################################################################################
########################                    BEGIN PIPELINE!                     ########################
########################################################################################################

# Check if input read files exist
if [ ! -s "$reads_1" ]; then echo "$reads_1 file does not exist. Exiting..."; exit 1; fi
if [ ! -s "$reads_2" ]; then echo "$reads_2 file does not exist. Exiting..."; exit 1; fi

if [ "$pre_qc_report" = true ]; then
    ########################################################################################################
    ########################                 MAKING PRE-QC REPORT                   ########################
    ########################################################################################################
    echo "Making Pre-QC report"
    mkdir -p "${out}/pre-QC_report"
    fastqc -q -t "$threads" -o "${out}/pre-QC_report" -f fastq "$reads_1" "$reads_2"

    if [ $? -ne 0 ]; then
        echo "Something went wrong with making pre-QC fastqc report. Exiting."
        exit 1
    fi
    rm "${out}/pre-QC_report"/*zip
    echo "Pre-QC report saved to: ${out}/pre-QC_report"
fi

########################################################################################################
########################                 RUNNING BBDuk STEPS                      ########################
########################################################################################################

# Step 1: Adapter trimming
echo "Running BBDuk step 1: Adapter trimming"
bbduk.sh -Xmx16g \
    in1="$reads_1" \
    in2="$reads_2" \
    out1="${out}/step1_adptrim_1.fastq" \
    out2="${out}/step1_adptrim_2.fastq" \
    ref="$ref" \
    ktrim=r k=23 mink=11 hdist=1 tpe tbo \
    minlen="$minlen" \
    threads="$threads" \
    gzip=false

if [ $? -ne 0 ]; then
    echo "Something went wrong with BBDuk step 1. Exiting."
    exit 1
fi

# Step 2: Quality trimming
echo "Running BBDuk step 2: Quality trimming"
bbduk.sh -Xmx16g \
    in1="${out}/step1_adptrim_1.fastq" \
    in2="${out}/step1_adptrim_2.fastq" \
    out1="${out}/step2_qualtrim_1.fastq" \
    out2="${out}/step2_qualtrim_2.fastq" \
    qtrim=r trimq="$trimq" \
    ftm=5 \
    minlen="$minlen" \
    threads="$threads" \
    gzip=false

if [ $? -ne 0 ]; then
    echo "Something went wrong with BBDuk step 2. Exiting."
    exit 1
fi

# Step 3: Trimming left bases
echo "Running BBDuk step 3: Trimming left $ftl bases"
bbduk.sh -Xmx16g \
    in1="${out}/step2_qualtrim_1.fastq" \
    in2="${out}/step2_qualtrim_2.fastq" \
    out1="${out}/step3_lefttrim_1.fastq" \
    out2="${out}/step3_lefttrim_2.fastq" \
    ftl="$ftl" \
    threads="$threads" \
    gzip=false

if [ $? -ne 0 ]; then
    echo "Something went wrong with BBDuk step 3. Exiting."
    exit 1
fi

# Optional deduplication with clumpify
if [ "$dedup" = true ]; then
    echo "Running Clumpify for deduplication"
    clumpify.sh \
        in1="${out}/step3_lefttrim_1.fastq" \
        in2="${out}/step3_lefttrim_2.fastq" \
        out1="${out}/final_reads_1.fastq" \
        out2="${out}/final_reads_2.fastq" \
        dedupe \
        threads="$threads" \
        gzip=false

    if [ $? -ne 0 ]; then
        echo "Something went wrong with Clumpify. Exiting."
        exit 1
    fi
else
    mv "${out}/step3_lefttrim_1.fastq" "${out}/final_reads_1.fastq"
    mv "${out}/step3_lefttrim_2.fastq" "${out}/final_reads_2.fastq"
fi

echo "Final processed reads are stored in: ${out}/final_reads_1.fastq and ${out}/final_reads_2.fastq"

if [ "$post_qc_report" = true ]; then
    ########################################################################################################
    ########################                 MAKING POST-QC REPORT                  ########################
    ########################################################################################################
    echo "Making Post-QC report"
    mkdir -p "${out}/post-QC_report"
    fastqc -t "$threads" -o "${out}/post-QC_report" -f fastq "${out}/final_reads_1.fastq" "${out}/final_reads_2.fastq"

    if [ $? -ne 0 ]; then
        echo "Something went wrong with making post-QC fastqc report. Exiting."
        exit 1
    fi
    rm "${out}/post-QC_report"/*zip
    echo "Post-QC report saved to: ${out}/post-QC_report"
fi

########################################################################################################
########################              READ QC PIPELINE COMPLETE!!!              ########################
########################################################################################################
echo "READ QC PIPELINE COMPLETE!!!"
