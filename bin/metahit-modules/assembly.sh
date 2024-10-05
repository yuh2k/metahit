"""
Megahit cannot be detected in the environment
metaspades can work
megahit parameters as options
metaspades with option -k 21,33,55,61 hint must be odd and smaller than 128
not use relative path
"""
#!/usr/bin/env bash

help_message () {
    echo ""
    echo "Usage: metahit assembly [options] -1 reads_1.fastq -2 reads_2.fastq -o output_dir"
    echo "Options:"
    echo "    -1 STR          forward fastq reads"
    echo "    -2 STR          reverse fastq reads"
    echo "    -o STR          output directory"
    echo "    -m INT          memory in GB (default=24)"
    echo "    -t INT          number of threads (default=1)"
    echo "    -l INT          minimum length of assembled contigs (default=1000)"
    echo ""
    echo "    --megahit       assemble with megahit (default)"
    echo "    --metaspades    assemble with metaspades instead of megahit"
    echo ""
}
SOFT="bin/metahit-scripts"

# default parameters
mem=24; threads=1; out="false"; reads_1="false"; reads_2="false"; min_len=1000
metaspades_assemble=false; megahit_assemble=true

# load in parameters
OPTS=$(getopt -o ht:m:o:1:2:l: --long help,metaspades,megahit -- "$@")
if [ $? != 0 ]; then help_message; exit 1; fi

eval set -- "$OPTS"

while true; do
    case "$1" in
        -t) threads=$2; shift 2;;
        -m) mem=$2; shift 2;;
        -o) out=$2; shift 2;;
        -1) reads_1=$2; shift 2;;
        -2) reads_2=$2; shift 2;;
        -l) min_len=$2; shift 2;;
        -h | --help) help_message; exit 1; shift 1;;
        --megahit) megahit_assemble=true; shift 1;;
        --metaspades) metaspades_assemble=true; shift 1;;
        --) shift; break ;;
        *) break;;
    esac
done

if [ "$out" = "false" ] || [ "$reads_1" = "false" ] || [ "$reads_2" = "false" ]; then 
    help_message; exit 1
fi

mkdir -p $out || { echo "Error: Cannot create output directory."; exit 1; }

if [ "$metaspades_assemble" = true ]; then
    echo "ASSEMBLING WITH METASPADES"
    mkdir -p ${out}/metaspades.tmp
    metaspades.py --tmp-dir ${out}/metaspades.tmp -k 21,33,55,61 -t $threads -m $mem -o ${out}/metaspades -1 $reads_1 -2 $reads_2
    if [ -d "${out}/metaspades.tmp" ]; then rm -r ${out}/metaspades.tmp; fi
    if [ ! -f "${out}/metaspades/scaffolds.fasta" ]; then echo "Error: metaSPAdes assembly failed."; exit 1; fi
fi

if [ "$megahit_assemble" = true ]; then
    echo "ASSEMBLING WITH MEGAHIT"
    megahit -1 $reads_1 -2 $reads_2 -o ${out}/megahit --min-contig-len 1000 --k-min 21 --k-max 141 --k-step 12 --merge-level 20,0.95 -t $threads -m ${mem}000000000
    if [ ! -f "${out}/megahit/final.contigs.fa" ]; then echo "Error: MEGAHIT assembly failed."; exit 1; fi
fi

if [ "$megahit_assemble" = true ] && [ "$metaspades_assemble" = true ]; then
    ${SOFT}/fix_megahit_contig_naming.py ${out}/megahit/final.contigs.fa $min_len > ${out}/megahit_final_assembly.fasta
    ${SOFT}/rm_short_contigs.py $min_len ${out}/metaspades/scaffolds.fasta > ${out}/metaspades_final_assembly.fasta
elif [ "$metaspades_assemble" = true ]; then
    ${SOFT}/rm_short_contigs.py $min_len ${out}/metaspades/scaffolds.fasta > ${out}/final_assembly.fasta
elif [ "$megahit_assemble" = true ]; then 
    ${SOFT}/fix_megahit_contig_naming.py ${out}/megahit/final.contigs.fa $min_len > ${out}/final_assembly.fasta
fi

if [[ ! -s ${out}/final_assembly.fasta ]]; then echo "Error: Final assembly failed."; exit 1; fi

echo "RUNNING ASSEMBLY QC WITH QUAST"
quast.py -t $threads -o ${out}/QUAST_out -m 500 ${out}/final_assembly.fasta
cp ${out}/QUAST_out/report.html ${out}/assembly_report.html

if [[ ! -s ${out}/assembly_report.html ]]; then echo "Error: QUAST analysis failed."; exit 1; fi

echo "ASSEMBLY PIPELINE COMPLETED SUCCESSFULLY!"
