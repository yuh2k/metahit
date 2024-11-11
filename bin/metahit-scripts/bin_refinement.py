#!/usr/bin/env python
import os
os.environ['NUMBA_DISABLE_JIT'] = '1'


import os
import argparse
import logging
import shutil
import scipy.sparse as scisp
import sys
from MetaCC.Script.normalized_contact import NormCCMap
from MetaCC.Script.cluster import ClusterBin
from MetaCC.Script.utils import gen_bins
from bin3C_python3.mzd.cluster import cluster_map, cluster_report, write_fasta
from ImputeCC.Script.imputation import ImputeMatrix
from ImputeCC.Script.pre_clustering import PreCluster
from ImputeCC.Script.final_cluster import FinalCluster
import pandas as pd

import argparse
import pandas as pd

def load_contig_info(file_path):
    print(f"Loading contig info from: {file_path}")
    try:
        contig_info = pd.read_csv(file_path)
        print(f"Successfully loaded contig info with shape: {contig_info.shape}")
        return contig_info
    except Exception as e:
        print(f"Failed to load contig_info.csv: {e}")
        raise

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run Binning Refinement")
    parser.add_argument('--method', choices=['metacc', 'bin3c', 'imputecc'], required=True, help='Choose refinement method')
    parser.add_argument('--contig_file', required=True, help='Path to contig_info.csv')
    parser.add_argument('--hic_matrix', required=True, help='Path to Hi-C normalized matrix (.npz)')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--fasta', help='Path to reference FASTA (required for bin3C)')
    parser.add_argument('--num_gene', type=int, help='Number of marker genes detected')
    parser.add_argument('--seed', type=int, default=None, help='Random seed')
    return parser.parse_args()

# Parse the arguments
args = parse_arguments()


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('refinement')

def run_metacc_refinement(contig_file, hic_matrix, outdir, num_gene=None, seed=None):
    logger.info("Starting MetaCC refinement...")
    
    # Step 1: Load contig_info.csv
    try:
        contig_info = pd.read_csv(contig_file, names=['name', 'sites', 'length', 'covcc'])
        if contig_info.empty or contig_info.isnull().values.any():
            raise ValueError("contig_info.csv is empty or contains NaN values")
        logger.info(f"Loaded contig_info with shape: {contig_info.shape}")
    except Exception as e:
        logger.error(f"Failed to load contig_info.csv: {e}")
        sys.exit(1)
    
    # Step 2: Load the Hi-C contact matrix from .npz file
    try:
        logger.info(f"Loading Hi-C matrix from: {hic_matrix}")
        seq_map = scisp.load_npz(hic_matrix)
        if seq_map is None or seq_map.shape[0] == 0:
            raise ValueError("Loaded contact matrix is empty or None")
        logger.info(f"Loaded contact matrix with shape: {seq_map.shape}")
    except Exception as e:
        logger.error(f"Failed to load Hi-C matrix: {e}")
        sys.exit(1)
    
    # Step 3: Initialize NormCCMap
    try:
        norm_result = None  # Replace this with actual norm_result calculation if needed
        thres = 5  # Example threshold value, adjust as needed
        norm_result_obj = NormCCMap(outdir, contig_info, seq_map, norm_result, thres)
    except Exception as e:
        logger.error(f"Failed to initialize NormCCMap: {e}")
        sys.exit(1)

    # Step 4: Perform clustering
    try:
        cluster_process = ClusterBin(outdir, norm_result_obj.name, norm_result_obj.len, norm_result_obj.seq_map, num_gene=num_gene, seed=seed)
        gen_bins(contig_file, os.path.join(outdir, 'cluster.txt'), os.path.join(outdir, 'BIN'))
    except Exception as e:
        logger.error(f"Clustering process failed: {e}")
        sys.exit(1)
    
    logger.info("MetaCC refinement finished.")


def run_bin3c_refinement(hic_matrix, outdir, fasta_file, seed=None):
    logger.info("Starting bin3C refinement...")
    clustering = cluster_map(hic_matrix, method='infomap', seed=seed, work_dir=outdir)
    cluster_report(hic_matrix, clustering)
    write_fasta(hic_matrix, outdir, clustering, source_fasta=fasta_file)
    logger.info("bin3C refinement finished.")

def run_imputecc_refinement(contig_info, hic_matrix, outdir, num_gene=None, seed=None):
    logger.info("Starting ImputeCC refinement...")
    imputer = ImputeMatrix(contig_info, hic_matrix)
    bins, bin_of_contigs = PreCluster(imputer.marker_contig_counts, imputer.marker_contigs, imputer.contig_markers)
    cluster_process = FinalCluster(imputer.contig_info, imputer.contig_local, 
                                   imputer.dict_contigRev, imputer.dict_contigRevLocal, 
                                   imputer.dict_contig_len, imputer.contig_markers, bins, bin_of_contigs,
                                   imputer.normcc_matrix, imputer.imputed_matrix)
    gen_bins(contig_info, os.path.join(outdir, 'cluster_imputecc.txt'), os.path.join(outdir, 'FINAL_BIN'))
    logger.info("ImputeCC refinement finished.")

def main():
    parser = argparse.ArgumentParser(description="Run Binning Refinement")
    parser.add_argument('--method', choices=['metacc', 'bin3c', 'imputecc'], required=True, help='Choose refinement method')
    parser.add_argument('--contig_file', required=True, help='Path to contig_info.csv')
    parser.add_argument('--hic_matrix', required=True, help='Path to Hi-C normalized matrix (.npz)')
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('--fasta', help='Path to reference FASTA (required for bin3C)')
    parser.add_argument('--num_gene', type=int, help='Number of marker genes detected')
    parser.add_argument('--seed', type=int, default=None, help='Random seed')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    if args.method == 'metacc':
        run_metacc_refinement(args.contig_file, args.hic_matrix, args.outdir, args.num_gene, args.seed)
    elif args.method == 'bin3c':
        if not args.fasta:
            raise ValueError("FASTA file is required for bin3C refinement")
        run_bin3c_refinement(args.hic_matrix, args.outdir, args.fasta, args.seed)
    elif args.method == 'imputecc':
        run_imputecc_refinement(args.contig_file, args.hic_matrix, args.outdir, args.num_gene, args.seed)

if __name__ == "__main__":
    # Ensure args is defined before using it
    args = parse_arguments()

    # Load the contig information
    contig_info = load_contig_info(args.contig_file)
    print("Loaded contig information successfully.")

