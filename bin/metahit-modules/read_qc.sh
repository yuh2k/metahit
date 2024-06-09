#!/usr/bin/env bash

###########################################################################################################################################################
#       	                                                                                                                                          #
# This script is meant to be a comprehensive solution to QC new HiSeq reads in preparation for assembly and other operations.                             #
# The main things this pipeline accomplishes are read trimming based on quality scores, and removal of human sequences.                                   #
# The script also produces a FASTQC report before and after the procedures.                                                                               #
#	                                                                                                                                 		  #
###########################################################################################################################################################

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
    echo "    -x STR          prefix of host index in bbmap database folder (default=hg38)"
    echo ""
    echo "    --skip-bbduk        don't remove human sequences with bbduk"
    echo "    --skip-trimming     don't trim sequences with bbduk"
    echo "    --skip-pre-qc-report don't make FastQC report of input sequences"
    echo "    --skip-post-qc-report don't make FastQC report of final sequences"
    echo "";}

########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################

# Default parameters
threads=1; out="false"; reads_1="false"; reads_2="false"
bbduk=true; trim=true; pre_qc_report=true; post_qc_report=true
HOST=hg38

# Load in parameters using getopt
OPTS=`getopt -o ht:o:1:2:x: --long help,skip-trimming,skip-bbduk,skip-pre-qc-report,skip-post-qc-report -- "$@"`
# Make sure the parameters are entered correctly
if [ $? -ne 0 ]; then help_message; exit 1; fi

# Loop through input parameters
while true; do
    case "$1" in
        -t) threads=$2; shift 2;;
        -o) out=$2; shift 2;;
        -1) reads_1=$2; shift 2;;
        -2) reads_2=$2; shift 2;;
        -x) HOST=$2; shift 2;;
        -h | --help) help_message; exit 1; shift 1;;
        --skip-trimming) trim=false; shift 1;;
        --skip-bbduk) bbduk=false; shift 1;;
        --skip-pre-qc-report) pre_qc_report=false; shift 1;;
        --skip-post-qc-report) post_qc_report=false; shift 1;;
        --) help_message; exit 1; shift; break ;;
        *) break;;
    esac
done

########################################################################################################
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################

# Check if all required parameters are entered
if [ "$out" = "false" ] || [ "$reads_1" = "false" ] || [ "$reads_2" = "false" ]; then 
    help_message; exit 1
fi

# Check if forward and reverse reads are the same file
if [ "$reads_1" = "$reads_2" ]; then
    echo "The forward and reverse reads are the same file. Exiting pipeline."
    exit 1
fi

# Check if output directory exists, if not, create it
if [ ! -d $out ]; then
        mkdir -p $out
else
        echo "Warning: $out already exists."
fi

########################################################################################################
########################                    BEGIN PIPELINE!                     ########################
########################################################################################################

# Check if input read files exist
if [ ! -s $reads_1 ]; then echo "$reads_1 file does not exist. Exiting..."; exit 1; fi
if [ ! -s $reads_2 ]; then echo "$reads_2 file does not exist. Exiting..."; exit 1; fi

if [ "$pre_qc_report" = true ]; then
    ########################################################################################################
    ########################                 MAKING PRE-QC REPORT                   ########################
    ########################################################################################################
    echo "Making Pre-QC report"
    mkdir -p ${out}/pre-QC_report
    fastqc -q -t $threads -o ${out}/pre-QC_report -f fastq $reads_1 $reads_2
    
    if [ $? -ne 0 ]; then 
        echo "Something went wrong with making pre-QC fastqc report. Exiting."
        exit 1
    fi
    rm ${out}/pre-QC_report/*zip
    echo "Pre-QC report saved to: ${out}/pre-QC_report"
fi

if [ "$trim" = true ]; then
    ########################################################################################################
    ########################                 RUNNING BBDuk FOR TRIMMING                    ########################
    ########################################################################################################
    echo "Running BBDuk for trimming"
    /Users/bach/Downloads/bbmap/bbduk.sh -Xmx16g in1=$reads_1 in2=$reads_2 out1=${out}/trimmed_1.fastq out2=${out}/trimmed_2.fastq ref=/Users/bach/Downloads/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=10 minlen=50
    
    if [ $? -ne 0 ]; then 
        echo "Something went wrong with BBDuk trimming. Exiting."
        exit 1
    fi
    echo "Trimmed reads saved to: ${out}/trimmed_1.fastq and ${out}/trimmed_2.fastq"
    
    reads_1=${out}/trimmed_1.fastq
    reads_2=${out}/trimmed_2.fastq
fi

if [ "$bbduk" = true ]; then
    ########################################################################################################
    ########################               REMOVING HOST SEQUENCES                  ########################
    ########################################################################################################
    echo "Removing host sequences with BBDuk"
    /Users/bach/Downloads/bbmap/bbduk.sh -Xmx16g in1=$reads_1 in2=$reads_2 out1=${out}/cleaned_1.fastq out2=${out}/cleaned_2.fastq ref=/Users/bach/Downloads/bbmap/resources/hg38.fa k=31 hdist=1
    
    if [ $? -ne 0 ]; then 
        echo "Something went wrong with BBDuk host removal. Exiting."
        exit 1
    fi
    echo "Cleaned reads saved to: ${out}/cleaned_1.fastq and ${out}/cleaned_2.fastq"
    
    reads_1=${out}/cleaned_1.fastq
    reads_2=${out}/cleaned_2.fastq
fi    

# Move final clean reads to their respective files
mv $reads_1 ${out}/final_pure_reads_1.fastq
mv $reads_2 ${out}/final_pure_reads_2.fastq
echo "Contamination-free and trimmed reads are stored in: ${out}/final_pure_reads_1.fastq and ${out}/final_pure_reads_2.fastq"

if [ "$post_qc_report" = true ]; then
    ########################################################################################################
    ########################                 MAKING POST-QC REPORT                  ########################
    ########################################################################################################
    echo "Making Post-QC report"
    mkdir -p ${out}/post-QC_report
    fastqc -t $threads -o ${out}/post-QC_report -f fastq ${out}/final_pure_reads_1.fastq ${out}/final_pure_reads_2.fastq
    
    if [ $? -ne 0 ]; then 
        echo "Something went wrong with making post-QC fastqc report. Exiting."
        exit 1
    fi
    rm ${out}/post-QC_report/*zip
    echo "Post-QC report saved to: ${out}/post-QC_report"
fi

########################################################################################################
########################              READ QC PIPELINE COMPLETE!!!              ########################
########################################################################################################
echo "READ QC PIPELINE COMPLETE!!!"