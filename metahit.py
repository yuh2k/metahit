#!/usr/bin/env python
import argparse
import subprocess
import os
import sys

script_dir=sys.path[0]

def run_command(command):
    try:
        print(f"[INFO] Executing command:\n{command}")
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Command failed: {e}")
        exit(1)

def ensure_dir_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def absolute_path(path):
    return os.path.abspath(path)

def find_fastq(filepath):
    """Detect if a file exists with either .fastq or .fastq.gz extension."""
    for ext in ["", ".gz"]:
        if os.path.exists(filepath + ext):
            return filepath + ext
    raise FileNotFoundError(f"{filepath} or {filepath}.gz not found.")

def readqc(args):
    print("[INFO] Running Read QC")
    output_dir = absolute_path(args.output)
    ensure_dir_exists(output_dir)
    r1_path = find_fastq(args.r1)
    r2_path = find_fastq(args.r2)
    command = script_dir+f"/bin/metahit-modules/read_qc.sh -p {script_dir} -1 {r1_path} -2 {r2_path} -o {output_dir} -t {args.threads}"
    if args.xmx:
        command += f" --xmx {args.xmx}"
    if args.ftm:
        command += f" --ftm {args.ftm}"
    run_command(command)

def assembly(args):
    print("[INFO] Running Assembly")
    output_dir = absolute_path(args.output)
    ensure_dir_exists(output_dir)
    r1_path = find_fastq(args.r1)
    r2_path = find_fastq(args.r2)
    command = script_dir+f"/bin/metahit-modules/assembly.sh -p {script_dir} -1 {r1_path} -2 {r2_path} -o {output_dir} -m {args.memory} -t {args.threads}"
    if args.min_len:
        command += f" -l {args.min_len}"
    if args.megahit:
        command += f" --megahit --k-min {args.k_min} --k-max {args.k_max} --k-step {args.k_step} --merge-level {args.merge_level}"
    elif args.metaspades:
        command += f" --metaspades --k-list {args.k_list}"

    run_command(command)

def alignment(args):
    print("[INFO] Running Alignment")
    output_dir = absolute_path(args.output)
    ensure_dir_exists(output_dir)
    r1_path = find_fastq(args.r1)
    r2_path = find_fastq(args.r2)
    command = script_dir+f"/bin/metahit-modules/alignment.sh -p {script_dir} -r {absolute_path(args.ref)} -1 {r1_path} -2 {r2_path} -o {output_dir} --threads {args.threads}"
    if args.samtools_filter:
        command += f" --samtools-filter '{args.samtools_filter}'"
    run_command(command)

def coverage_estimation(args):
    print("[INFO] Running Coverage Estimation")
    output_dir = absolute_path(args.output)
    ensure_dir_exists(output_dir)
    r1_path = find_fastq(args.r1)
    r2_path = find_fastq(args.r2)
    command = script_dir+f"/bin/metahit-modules/coverage_estimation.sh -p {script_dir} -1 {r1_path} -2 {r2_path} -r {absolute_path(args.ref)} -o {output_dir}"
    run_command(command)

def raw_contact(args):
    print("[INFO] Running Raw Contact Generation")
    output_dir = absolute_path(args.output)
    ensure_dir_exists(output_dir)
    
    # Check if enzyme argument is provided
    if not args.enzyme:
        print("[ERROR] Missing enzyme argument.")
        exit(1)

    command = (
        script_dir+f"/bin/metahit-modules/raw_contact1.sh -p {script_dir} --bam {absolute_path(args.bam)} "
        f"--fasta {absolute_path(args.fasta)} --out {output_dir} --enzyme {args.enzyme}"
    )
    
    # Optional coverage file
    if args.coverage:
        command += f" --coverage {absolute_path(args.coverage)}"
    
    run_command(command)


def normalization(args):
    print(f"[INFO] Running {args.method} Normalization")
    output_dir = absolute_path(args.output)
    ensure_dir_exists(output_dir)

    # Construct the normalization command
    command = f'"{script_dir}/bin/metahit-modules/normalization.sh" {args.method} -p "{script_dir}" --contig_file "{absolute_path(args.contig_file)}" --contact_matrix_file "{absolute_path(args.contact_matrix_file)}" --output "{output_dir}" --thres {args.thres}'

    if args.method == "raw":
        command += f" --min_len {args.min_len} --min_signal {args.min_signal}"
    elif args.method == "normcc":
        print("[INFO] Running NormCC Normalization")
    elif args.method == "hiczin":
        command += f" --epsilon {args.epsilon}"
    elif args.method == "bin3c":
        command += f" --max_iter {args.max_iter} --tol {args.tol}"
    elif args.method == "metator":
        command += f" --epsilon {args.epsilon}"
    elif args.method == "fastnorm":
        command += f" --epsilon {args.epsilon}"

    print(f"[DEBUG] Running command: {command}")
    run_command(command)

def bin_refinement(args):
    print("[INFO] Running Bin Refinement")
    output_dir = absolute_path(args.output)
    ensure_dir_exists(output_dir)
    
    fasta_file = absolute_path(args.fasta)
    bam_file = absolute_path(args.bam)
    
    # Construct the bin_refinement.sh command
    command = script_dir+f"/bin/metahit-modules/bin_refinement.sh {fasta_file} {bam_file} {output_dir} {script_dir}"
    
    # Append optional arguments
    if args.threads:
        command += f" -t {args.threads}"
    if args.enzyme:
        for enzyme in args.enzyme:
            command += f" --enzyme {enzyme}"
    if args.metacc_min_len:
        command += f" --metacc-min-len {args.metacc_min_len}"
    if args.metacc_min_signal:
        command += f" --metacc-min-signal {args.metacc_min_signal}"
    if args.bin3c_min_len:
        command += f" --bin3c-min-len {args.bin3c_min_len}"
    if args.bin3c_min_signal:
        command += f" --bin3c-min-signal {args.bin3c_min_signal}"
    if args.thres:
        command += f" --thres {args.thres}"
    if args.cover:
        command += f" --cover"
    
    # Execute the command
    run_command(command)


def scaffolding(args):
    print("[INFO] Running Scaffolding")
    output_dir = os.path.abspath(args.outdir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    command = (
        f'"{script_dir}/bin/metahit-modules/scaffolding.sh" '
        f'-p "{script_dir}" '
        f'--fasta "{args.fasta}" --bam "{args.bam}" --enzyme "{args.enzyme}" '
        f'--hic1 "{args.hic1}" --hic2 "{args.hic2}" --outdir "{output_dir}" '
        f'-t {args.threads} -m {args.memory} -r {args.resolution}'
    )
    run_command(command)


def reassembly(args):
    print("[INFO] Running Reassembly")
    output_dir = os.path.abspath(args.outdir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    command = (
        f'"{script_dir}/bin/metahit-modules/reassembly.sh" '
        f'-p "{script_dir}" '
        f'--bin "{args.bin}" --hic1 "{args.hic1}" --hic2 "{args.hic2}" '
        f'--sg1 "{args.sg1}" --sg2 "{args.sg2}" --bam "{args.bam}" '
        f'--outdir "{output_dir}" -t {args.threads} -m {args.memory}'
    )
    run_command(command)


def viralcc(args):
    print("[INFO] Running ViralCC pipeline")
    output_dir = os.path.abspath(args.OUTDIR)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    command = (
        f'"{script_dir}/bin/metahit-modules/viral_binning.sh" pipeline '
        f'-p "{script_dir}" '
    )
    if args.min_len:
        command += f"--min-len {args.min_len} "
    if args.min_mapq:
        command += f"--min-mapq {args.min_mapq} "
    if args.min_match:
        command += f"--min-match {args.min_match} "
    if args.min_k:
        command += f"--min-k {args.min_k} "
    if args.random_seed:
        command += f"--random-seed {args.random_seed} "

    # 必选参数：FASTA, BAM, VIRAL, OUTDIR
    command += f'"{args.FASTA}" "{args.BAM}" "{args.VIRAL}" "{output_dir}"'

    run_command(command)


def annotation(args):
    print("[INFO] Running Annotation via GTDB-Tk classify_wf")
    output_dir = os.path.abspath(args.outdir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    command = (
        f'"{script_dir}/bin/metahit-modules/annotation.sh" '
        f'-p "{script_dir}" '
        f'--genome_dir "{args.genome_dir}" '
        f'--out_dir "{output_dir}" '
    )
    if args.extension:
        command += f'--extension {args.extension} '
    if args.cpus:
        command += f'--cpus {args.cpus} '

    run_command(command)



def main():
    parser = argparse.ArgumentParser(description="MetaHit Pipeline Command Line Interface")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Read QC
    qc_parser = subparsers.add_parser("readqc")
    qc_parser.add_argument("-1", dest="r1", required=True, help="Path to R1 reads file")
    qc_parser.add_argument("-2", dest="r2", required=True, help="Path to R2 reads file")
    qc_parser.add_argument("-o", "--output", required=True, help="Output directory")
    qc_parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads")
    qc_parser.add_argument("--xmx", help="Memory for Java in BBDuk (default is 60% of system memory)")
    qc_parser.add_argument("--ftm", type=int, help="ftm value for BBDuk (default=5)")
    qc_parser.add_argument("--minlen", type=int, default=50, help="Minimum read length after trimming (default=50)")
    qc_parser.add_argument("--trimq", type=int, default=10, help="Quality threshold for trimming (default=10)")
    qc_parser.add_argument("--ftl", type=int, default=10, help="Trim bases from the left (default=10)")
    qc_parser.add_argument("--dedup", action="store_true", help="Perform deduplication with Clumpify (default=False)")
    qc_parser.add_argument("--skip-pre-qc", action="store_true", help="Skip FastQC for input reads")
    qc_parser.add_argument("--skip-post-qc", action="store_true", help="Skip FastQC for final reads")
    qc_parser.add_argument("--k", type=int, help="k-mer size for adapter trimming (default=23)")
    qc_parser.add_argument("--mink", type=int, help="Minimum k-mer size (default=11)")
    qc_parser.add_argument("--hdist", type=int, help="Hamming distance for k-mer matching (default=1)")

    # Assembly
    asm_parser = subparsers.add_parser("assembly")
    asm_parser.add_argument("-1", dest="r1", required=True, help="Path to R1 reads file")
    asm_parser.add_argument("-2", dest="r2", required=True, help="Path to R2 reads file")
    asm_parser.add_argument("-o", "--output", required=True, help="Output directory")
    asm_parser.add_argument("-m", "--memory", type=int, default=24, help="Memory in GB")
    asm_parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads")
    asm_parser.add_argument("--megahit", action="store_true", help="Use MEGAHIT for assembly")
    asm_parser.add_argument("--metaspades", action="store_true", help="Use metaSPAdes for assembly")
    asm_parser.add_argument("--k-min", type=int, default=21, help="Minimum k-mer size")
    asm_parser.add_argument("--k-max", type=int, default=141, help="Maximum k-mer size")
    asm_parser.add_argument("--k-step", type=int, default=12, help="k-mer step size")
    asm_parser.add_argument("--k-list", default="21,33,55,77", help="List of k-mers for metaSPAdes")
    asm_parser.add_argument("-l", "--min-len", type=int, default=1000, help="Minimum contig length")
    asm_parser.add_argument("--merge-level", default="20,0.95", help="Merge level for MEGAHIT (default='20,0.95')")


    # Alignment
    aln_parser = subparsers.add_parser("alignment")
    aln_parser.add_argument("-r", "--ref", required=True, help="Path to reference assembly file")
    aln_parser.add_argument("-1", dest="r1", required=True, help="Path to R1 reads file")
    aln_parser.add_argument("-2", dest="r2", required=True, help="Path to R2 reads file")
    aln_parser.add_argument("-o", "--output", required=True, help="Output directory")
    aln_parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    aln_parser.add_argument("--samtools-filter", help="Samtools filter options")

    # Coverage Estimation
    cov_parser = subparsers.add_parser("coverage_estimation")
    cov_parser.add_argument("-1", dest="r1", required=True, help="Path to R1 reads file")
    cov_parser.add_argument("-2", dest="r2", required=True, help="Path to R2 reads file")
    cov_parser.add_argument("-r", "--ref", required=True, help="Path to reference assembly file")
    cov_parser.add_argument("-o", "--output", required=True, help="Output directory")

    # Raw Contact Generation
    raw_parser = subparsers.add_parser("raw_contact")
    raw_parser.add_argument("--bam", required=True, help="Path to BAM file")
    raw_parser.add_argument("--fasta", required=True, help="Path to FASTA file")
    raw_parser.add_argument("--out", dest="output", required=True, help="Output directory")
    raw_parser.add_argument("--coverage", help="Path to coverage file")
    raw_parser.add_argument("--enzyme", required=True, help="Enzyme name used in Hi-C experiment")


    # Normalization Command
    norm_parser = subparsers.add_parser("normalization")
    norm_parser.add_argument("method", choices=["raw", "normcc", "hiczin", "bin3c", "metator", "fastnorm"], help="Normalization method")
    norm_parser.add_argument("--contig_file", required=True, help="Path to contig info file")
    norm_parser.add_argument("--contact_matrix_file", required=True, help="Path to contact matrix file")
    norm_parser.add_argument("--output", required=True, help="Output directory")
    norm_parser.add_argument("--thres", type=int, default=5, help="Threshold value")
    norm_parser.add_argument("--min_len", type=int, default=500, help="Minimum length (raw normalization)")
    norm_parser.add_argument("--min_signal", type=int, default=1, help="Minimum signal (raw normalization)")
    norm_parser.add_argument("--max_iter", type=int, default=1000, help="Maximum iterations (bin3c normalization)")
    norm_parser.add_argument("--tol", type=float, default=1e-6, help="Tolerance (bin3c normalization)")
    norm_parser.add_argument("--epsilon", type=float, default=1, help="Epsilon value (fastnorm)")

    # Binning Refinement Command
    refinement_parser = subparsers.add_parser("bin_refinement", help="Perform bin refinement using bin_refinement.sh")
    refinement_parser.add_argument("--fasta", required=True, help="Path to the reference fasta sequence")
    refinement_parser.add_argument("--bam", required=True, help="Path to the input BAM file in query order")
    refinement_parser.add_argument("-o", "--output", required=True, help="Output directory for refined bins")
    refinement_parser.add_argument("-t", "--threads", type=int, default=10, help="Number of threads (default: 10)")
    refinement_parser.add_argument("-e", "--enzyme", action='append', help="Case-sensitive enzyme name. Use multiple times for multiple enzymes")
    refinement_parser.add_argument("--metacc-min-len", type=int, default=1000, help="Minimum acceptable contig length for metacc (default: 1000)")
    refinement_parser.add_argument("--metacc-min-signal", type=int, default=2, help="Minimum acceptable Hi-C signal for metacc (default: 2)")
    refinement_parser.add_argument("--bin3c-min-len", type=int, default=1000, help="Minimum acceptable contig length for bin3c (default: 1000)")
    refinement_parser.add_argument("--bin3c-min-signal", type=int, default=1, help="Minimum acceptable Hi-C signal for bin3c (default: 1)")
    refinement_parser.add_argument("--thres", type=float, default=0.01, help="Fraction threshold for NormCC-normalized Hi-C contacts (default: 0.01)")
    refinement_parser.add_argument("--cover", action='store_true', help="Overwrite existing files if set")

    # scaffolding子命令
    scaff_parser = subparsers.add_parser("scaffolding", help="Perform scaffolding")
    scaff_parser.add_argument("--fasta", required=True, help="Reference FASTA")
    scaff_parser.add_argument("--bam", required=True, help="Hi-C BAM file")
    scaff_parser.add_argument("--enzyme", required=True, help="Restriction enzyme")
    scaff_parser.add_argument("--hic1", required=True, help="Hi-C library forward fastq")
    scaff_parser.add_argument("--hic2", required=True, help="Hi-C library reverse fastq")
    scaff_parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    scaff_parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads")
    scaff_parser.add_argument("-m", "--memory", type=int, default=24, help="Memory in GB")
    scaff_parser.add_argument("-r", "--resolution", type=int, default=10000, help="Resolution (default 10kb)")

    # reassembly子命令
    reasm_parser = subparsers.add_parser("reassembly", help="Perform reassembly")
    reasm_parser.add_argument("--bin", required=True, help="Binning result directory")
    reasm_parser.add_argument("--hic1", required=True, help="Hi-C library forward fastq")
    reasm_parser.add_argument("--hic2", required=True, help="Hi-C library reverse fastq")
    reasm_parser.add_argument("--sg1", required=True, help="Shotgun forward fastq")
    reasm_parser.add_argument("--sg2", required=True, help="Shotgun reverse fastq")
    reasm_parser.add_argument("--bam", required=True, help="Hi-C BAM mapped to shotgun assembly")
    reasm_parser.add_argument("--outdir", required=True, help="Output directory")
    reasm_parser.add_argument("-p", "--metahit_path", required=False, help="Path to metahit folder", default=script_dir)
    reasm_parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads")
    reasm_parser.add_argument("-m", "--memory", type=int, default=24, help="Memory in GB")

    # ViralCC
    viralcc_parser = subparsers.add_parser("viralcc", help="Run ViralCC pipeline")
    viralcc_parser.add_argument("viralcc_subcommand", choices=["pipeline"], help="ViralCC command")
    viralcc_parser.add_argument("--min-len", type=int, help="Minimum acceptable reference length")
    viralcc_parser.add_argument("--min-mapq", type=int, help="Minimum acceptable mapping quality")
    viralcc_parser.add_argument("--min-match", type=int, help="Minimum acceptable matches")
    viralcc_parser.add_argument("--min-k", type=int, help="Lower bound of k")
    viralcc_parser.add_argument("--random-seed", type=int, help="Random seed")
    viralcc_parser.add_argument("FASTA", help="Reference fasta")
    viralcc_parser.add_argument("BAM", help="BAM file")
    viralcc_parser.add_argument("VIRAL", help="Viral contig labels")
    viralcc_parser.add_argument("OUTDIR", help="Output directory")

    # Annotation
    annotation_parser = subparsers.add_parser("annotation", help="Run annotation using GTDB-Tk classify_wf")
    annotation_parser.add_argument("--genome_dir", required=True, help="Directory containing bin genomes")
    annotation_parser.add_argument("--out_dir", dest="outdir", required=True, help="Output directory for annotation results")
    annotation_parser.add_argument("--extension", help="Genome file extension (default: fa)")
    annotation_parser.add_argument("--cpus", type=int, default=4, help="Number of CPUs (default: 4)")
    annotation_parser.set_defaults(func=annotation)


    # Link the subcommand to the function
    refinement_parser.set_defaults(func=bin_refinement)

    args = parser.parse_args()

    if args.command == "readqc":
        readqc(args)
    elif args.command == "assembly":
        assembly(args)
    elif args.command == "alignment":
        alignment(args)
    elif args.command == "coverage_estimation":
        coverage_estimation(args)
    elif args.command == "raw_contact":
        raw_contact(args)
    elif args.command == "normalization":
        normalization(args)
    elif args.command == "bin_refinement":
        bin_refinement(args)
    elif args.command == "scaffolding":
        scaffolding(args)
    elif args.command == "reassembly":
        reassembly(args)
    elif args.command == "viralcc":
        if args.viralcc_subcommand == "pipeline":
            viralcc(args)
    elif args.command == "annoation":
        annotation(args)



if __name__ == "__main__":
    main()
