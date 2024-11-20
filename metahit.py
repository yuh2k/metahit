#!/usr/bin/env python
import argparse
import subprocess
import os

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
    command = f"./bin/metahit-modules/read_qc.sh -1 {r1_path} -2 {r2_path} -o {output_dir} -t {args.threads}"
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
    command = f"./bin/metahit-modules/assembly.sh -1 {r1_path} -2 {r2_path} -o {output_dir} -m {args.memory} -t {args.threads}"
    if args.megahit:
        command += f" --megahit --k-min {args.k_min} --k-max {args.k_max} --k-step {args.k_step}"
    elif args.metaspades:
        command += f" --metaspades --k-list {args.k_list}"
    command += f" -l {args.min_len}"
    run_command(command)

def alignment(args):
    print("[INFO] Running Alignment")
    output_dir = absolute_path(args.output)
    ensure_dir_exists(output_dir)
    r1_path = find_fastq(args.r1)
    r2_path = find_fastq(args.r2)
    command = f"./bin/metahit-modules/alignment.sh -r {absolute_path(args.ref)} -1 {r1_path} -2 {r2_path} -o {output_dir} --threads {args.threads}"
    if args.samtools_filter:
        command += f" --samtools-filter '{args.samtools_filter}'"
    run_command(command)

def coverage_estimation(args):
    print("[INFO] Running Coverage Estimation")
    output_dir = absolute_path(args.output)
    ensure_dir_exists(output_dir)
    r1_path = find_fastq(args.r1)
    r2_path = find_fastq(args.r2)
    command = f"./bin/metahit-modules/coverage_estimation.sh -1 {r1_path} -2 {r2_path} -r {absolute_path(args.ref)} -o {output_dir}"
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
        f"./bin/metahit-modules/raw_contact1.sh --bam {absolute_path(args.bam)} "
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
    command = f"./bin/metahit-modules/normalization.sh {args.method} --contig_file {absolute_path(args.contig_file)} --contact_matrix_file {absolute_path(args.contact_matrix_file)} --output {output_dir} --thres {args.thres}"

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
    command = f"./bin/metahit-modules/bin_refinement.sh {fasta_file} {bam_file} {output_dir}"
    
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

if __name__ == "__main__":
    main()
