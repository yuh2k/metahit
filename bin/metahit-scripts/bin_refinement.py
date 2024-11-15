#!/usr/bin/env python3

import os
import argparse
import logging
from MetaCC.Script.utils import load_object, save_object, make_dir, gen_bins
from MetaCC.Script.cluster import ClusterBin
from ImputeCC.Script.imputation import ImputeMatrix
from ImputeCC.Script.final_cluster import FinalCluster
from ImputeCC.Script.pre_clustering import PreCluster
from bin3C_python3.mzd.cluster import cluster_map, write_report, write_fasta
from Bio import SeqIO

def setup_logging(output_dir):
    log_path = os.path.join(output_dir, 'bin_refinement.log')
    logging.basicConfig(filename=log_path, level=logging.DEBUG,
                        format='%(levelname)s: %(message)s')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    logging.getLogger('').addHandler(console)
    return logging.getLogger()

def generate_bins(args, logger):
    # 加载 contig 信息和 Hi-C 矩阵
    contig_info = load_object(args.contig_file)
    hic_matrix = load_object(args.contact_matrix_file)
    
    # Step 1: 生成 MetaCC Bins
    logger.info('Generating bins using MetaCC...')
    cluster_process = ClusterBin(args.output, contig_info['name'], contig_info['len'],
                                 hic_matrix, args.min_binsize, args.num_gene, args.seed)
    gen_bins(args.fasta, os.path.join(args.output, 'metacc_bins'), os.path.join(args.output, 'BIN'))

    # Step 2: 使用 Bin3C 方法生成 Bins
    logger.info('Generating bins using Bin3C...')
    clustering = cluster_map(hic_matrix, method='infomap', seed=args.seed, work_dir=args.output)
    write_report(hic_matrix, clustering, is_spades=not args.no_spades)
    write_fasta(hic_matrix, args.output, clustering, source_fasta=args.fasta, clobber=True)

    # Step 3: 使用 ImputeCC 方法生成 Bins
    logger.info('Generating bins using ImputeCC...')
    impute_matrix = ImputeMatrix(contig_info, hic_matrix, args.marker_file,
                                 gene_cov=args.gene_cov, rwr_rp=args.rwr_rp, rwr_thres=args.rwr_thres)
    save_object(os.path.join(args.output, 'ImputeCC_storage'), impute_matrix)
    pre_bins, _ = PreCluster(impute_matrix.marker_contig_counts, impute_matrix.marker_contigs,
                             impute_matrix.contig_markers, impute_matrix.imputed_matrix)
    final_cluster = FinalCluster(impute_matrix.contig_info, pre_bins)
    gen_bins(args.fasta, os.path.join(args.output, 'imputecc_bins'), os.path.join(args.output, 'FINAL_BIN'))

def refine_bins(args, logger):
    # 使用 MetaWRAP 进行 Bin Refinement
    logger.info('Refining bins using MetaWRAP...')
    metawrap_cmd = f"metaWRAP bin_refinement -t {args.threads} -c 50 -o {args.output}/refined_bins " \
                   f"-A {args.output}/BIN -B {args.output}/FINAL_BIN -C {args.output}/fasta"
    os.system(metawrap_cmd)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', required=True, help='Input FASTA file')
    parser.add_argument('--contig_file', required=True, help='Path to contig info file')
    parser.add_argument('--contact_matrix_file', required=True, help='Path to Hi-C contact matrix file')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--marker_file', help='Path to marker gene file')
    parser.add_argument('--min_binsize', type=int, default=150000, help='Minimum bin size')
    parser.add_argument('--num_gene', type=int, help='Number of marker genes')
    parser.add_argument('--threads', type=int, default=30, help='Number of threads')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    parser.add_argument('--gene_cov', type=float, default=0.9, help='Gene coverage')
    parser.add_argument('--rwr_rp', type=float, default=0.5, help='Random walk restart probability')
    parser.add_argument('--rwr_thres', type=int, default=80, help='RWR threshold')
    parser.add_argument('--no_spades', action='store_true', help='Assembly was not done using SPAdes')

    args = parser.parse_args()

    # Setup logging
    logger = setup_logging(args.output)

    try:
        generate_bins(args, logger)
        refine_bins(args, logger)
        logger.info('Bin generation and refinement completed successfully.')
    except Exception as e:
        logger.error(f'Error during bin refinement: {str(e)}')

if __name__ == '__main__':
    main()