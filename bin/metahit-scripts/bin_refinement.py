#!/usr/bin/env python
import os
import sys
import logging
import shutil
import argparse
import warnings
import subprocess

# Import necessary modules for each binning method
from MetaCC.Script.cluster import ClusterBin
from MetaCC.Script.predict_species_number import gen_bestk
from MetaCC.Script.utils import gen_bins, make_random_seed, make_dir
from bin3C_python3.mzd.cluster import cluster_map, cluster_report, write_fasta, write_report
from ImputeCC.Script.imputation import ImputeMatrix
from ImputeCC.Script.final_cluster import FinalCluster
from ImputeCC.Script.pre_clustering import PreCluster
from ImputeCC.Script.pre_common import get_bin_dirs
from ImputeCC.Script.pre_profile import Profile
from ImputeCC.Script.utility import save_object, load_object
from MetaCC.Script.normalized_contact import NormCCMap
from MetaCC.Script.normcc import normcc

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('binning')

# Default parameters for each method
metacc_runtime_defaults = {
    'min_binsize': 150000,
    'thres': 0.05,
}

imputecc_runtime_defaults = {
    'gene_cov': 0.9,
    'rwr_rp': 0.5,
    'rwr_thres': 80,
    'max_markers': 8000,
    'intra': 50,
    'inter': 0,
    'min_binsize': 100000,
    'cont_weight': 2,
    'min_comp': 50.0,
    'max_cont': 10.0,
    'report_quality': 10.0
}

def ifelse(arg, default):
    return default if arg is None else arg

def replace_in_file(file_path):
    """Function to adjust the naming convention for bin3C output FASTA files."""
    with open(file_path, 'r') as file:
        lines = file.readlines()
    with open(file_path, 'w') as file:
        for line in lines:
            if ">" in line and ":" in line:
                new_line = line.split(' ')[1].split(':')[1]
                file.write(">" + new_line + "\n")
            else:
                file.write(line)

def run_contact_map_generation(bam_file, enzyme_list, fasta_file, output_dir, metacc_folder, bin3c_folder):
    """Generate contact maps for MetaCC and bin3C."""
    from scripts.raw_contact_both import ContactMatrix
    logger.info('Generating contact maps for MetaCC and bin3C...')

    # Create necessary directories
    make_dir(metacc_folder)
    make_dir(bin3c_folder)

    # Initialize ContactMatrix object
    cm = ContactMatrix(
        bam_file,
        enzyme_list,
        fasta_file,
        output_dir,
        metacc_folder,
        bin3c_folder,
        min_insert=0,
        bin_size=None,
        # Parameters for MetaCC
        min_mapq_metacc=30,
        min_len_metacc=1000,
        min_match_metacc=30,
        min_signal_metacc=2,
        # Parameters for bin3C
        min_mapq_bin3c=60,
        min_len_bin3c=1000,
        min_match_bin3c=10,
        min_signal_bin3c=5
    )

    logger.info('Contact maps generated.')
    # Save the ContactMatrix object
    save_object(os.path.join(metacc_folder, 'contact_map.p'), cm)
    save_object(os.path.join(bin3c_folder, 'contact_map.p'), cm)

def run_normcc_normalization(metacc_folder):
    """Perform NormCC normalization."""
    logger.info('Performing NormCC normalization for MetaCC...')
    # Paths
    contig_file = os.path.join(metacc_folder, 'contig_info.csv')
    contact_matrix_file = os.path.join(metacc_folder, 'contact_matrix_metacc.npz')
    output_path = metacc_folder

    # Run NormCC normalization
    norm_result = normcc(contig_file, contact_matrix_file, output_path, thres=metacc_runtime_defaults['thres'])

    # Save the NormCCMap object
    hzmap = NormCCMap(
        seq_info=None,  # Will be loaded from files
        seq_map=None,
        norm_result=norm_result
    )
    save_object(os.path.join(metacc_folder, 'NormCC_normalized_contact.gz'), hzmap)

    logger.info('NormCC normalization completed.')

def run_metacc_binning(metacc_folder, fasta_file, num_gene, seed):
    """MetaCC binning process."""
    logger.info('Starting MetaCC binning process...')
    
    # Load the NormCCMap object
    normcc_object_path = os.path.join(metacc_folder, 'NormCC_normalized_contact.gz')
    logger.info(f'Loading NormCCMap object from {normcc_object_path}')
    
    hzmap = load_object(normcc_object_path)

    # Verify that the object loaded is actually of type NormCCMap
    if not isinstance(hzmap, NormCCMap):
        logger.error(f"Loaded object type: {type(hzmap)}")
        raise TypeError("Loaded object is not of type 'NormCCMap'")

    # If num_gene is not provided, generate it
    if num_gene is None:
        logger.info('Scanning marker genes...')
        num_gene = gen_bestk(metacc_folder, fasta_file)
        if num_gene == 0:
            logger.warning('No marker genes detected in the assembled contigs.')
    
    logger.info(f'{num_gene} marker genes detected.')

    if not seed:
        seed = make_random_seed()
    
    logger.info(f'Using random seed: {seed}')
    
    # Perform clustering
    cluster_process = ClusterBin(
        metacc_folder,
        hzmap.name,
        hzmap.len,
        hzmap.seq_map,
        metacc_runtime_defaults['min_binsize'],
        num_gene,
        seed
    )
    
    logger.info('Writing MetaCC bins...')
    # Generate bins
    bin_output_dir = os.path.join(metacc_folder, 'BIN')
    os.makedirs(bin_output_dir, exist_ok=True)
    cluster_file = os.path.join(metacc_folder, 'tmp', 'cluster.txt')
    gen_bins(fasta_file, cluster_file, bin_output_dir)
    
    logger.info('MetaCC binning completed.')

def run_bin3c_binning(bin3c_folder, fasta_file, seed):
    """bin3C binning process."""
    logger.info('Starting bin3C binning process...')
    
    # Load the ContactMatrix object
    cm_path = os.path.join(bin3c_folder, 'contact_map.p')
    logger.info(f'Loading ContactMatrix object from {cm_path}')
    cm = load_object(cm_path)

    clustering = cluster_map(cm, method='infomap', seed=seed, work_dir=bin3c_folder)
    cluster_report(cm, clustering)
    save_object(os.path.join(bin3c_folder, 'clustering.p'), clustering)
    
    write_report(os.path.join(bin3c_folder, 'cluster_report.csv'), clustering)

    write_fasta(cm, bin3c_folder, clustering, source_fasta=fasta_file, clobber=True, only_large=False)
    
    # Adjust the contig naming in the output FASTA files
    fasta_dir = os.path.join(bin3c_folder, 'fasta')
    for filename in os.listdir(fasta_dir):
        if filename.endswith('.fna'):
            file_path = os.path.join(fasta_dir, filename)
            replace_in_file(file_path)
    
    logger.info('bin3C binning completed.')

def run_imputecc_binning(imputecc_folder, contig_info, hic_matrix, fasta_file, num_gene):
    """ImputeCC binning process."""
    logger.info('Starting ImputeCC binning process...')
    
    marker_file = os.path.join(imputecc_folder, 'contigs.hmmout')
    if not os.path.exists(marker_file):
        num_gene = gen_bestk(imputecc_folder, fasta_file)
    
    # Create ImputeMatrix object
    imp = ImputeMatrix(
        contig_info,
        hic_matrix,
        marker_file,
        gene_cov=imputecc_runtime_defaults['gene_cov'],
        rwr_rp=imputecc_runtime_defaults['rwr_rp'],
        rwr_thres=imputecc_runtime_defaults['rwr_thres'],
        max_markers=imputecc_runtime_defaults['max_markers']
    )
    # Save the ImputeMatrix object
    save_object(os.path.join(imputecc_folder, 'ImputeCC_storage'), imp)
    
    bins, bin_of_contigs = PreCluster(
        imp.marker_contig_counts,
        imp.marker_contigs,
        imp.contig_markers,
        imp.imputed_matrix,
        imp.dict_contigRevLocal,
        intra=imputecc_runtime_defaults['intra'],
        inter=imputecc_runtime_defaults['inter']
    )
    
    # Final clustering
    checkm_bac_gene_table = os.path.join(imputecc_folder, 'checkm_bac_gene_table.tsv')
    checkm_ar_gene_table = os.path.join(imputecc_folder, 'checkm_ar_gene_table.tsv')
    
    cluster_process = FinalCluster(
        imp.contig_info,
        imp.contig_local,
        imp.dict_contigRev,
        imp.dict_contigRevLocal,
        imp.dict_contig_len,
        imp.contig_markers,
        bins,
        bin_of_contigs,
        imp.normcc_matrix,
        imp.imputed_matrix,
        checkm_bac_gene_table,
        checkm_ar_gene_table,
        imputecc_folder,
        intra=imputecc_runtime_defaults['intra'],
        inter=imputecc_runtime_defaults['inter'],
        cont_weight=imputecc_runtime_defaults['cont_weight'],
        min_comp=imputecc_runtime_defaults['min_comp'],
        max_cont=imputecc_runtime_defaults['max_cont'],
        report_quality=imputecc_runtime_defaults['report_quality'],
        min_binsize=imputecc_runtime_defaults['min_binsize']
    )
    
    logger.info('Writing ImputeCC bins...')
    bin_output_dir = os.path.join(imputecc_folder, 'FINAL_BIN')
    os.makedirs(bin_output_dir, exist_ok=True)
    cluster_file = os.path.join(imputecc_folder, 'cluster_imputecc.txt')
    gen_bins(fasta_file, cluster_file, bin_output_dir)
    logger.info('ImputeCC binning completed.')

def main():
    parser = argparse.ArgumentParser(description="Run Binning Processes")
    parser.add_argument('--bam', required=True, help='Path to BAM file')
    parser.add_argument('--enzyme', required=True, nargs='+', help='Enzyme list')
    parser.add_argument('--fasta', required=True, help='Path to reference FASTA')
    parser.add_argument('--coverage', help='Path to coverage file')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--contig_file', help='Path to contig_info.csv (required for ImputeCC)')
    parser.add_argument('--num_gene', type=int, help='Number of marker genes detected')
    parser.add_argument('--seed', type=int, default=None, help='Random seed')
    args = parser.parse_args()

    # Output folders
    metacc_folder = os.path.join(args.output, 'metacc')
    bin3c_folder = os.path.join(args.output, 'bin3c')
    imputecc_folder = os.path.join(args.output, 'imputecc')

    # Step 1: Generate Contact Maps
    run_contact_map_generation(
        bam_file=args.bam,
        enzyme_list=args.enzyme,
        fasta_file=args.fasta,
        output_dir=args.output,
        metacc_folder=metacc_folder,
        bin3c_folder=bin3c_folder
    )

    # Step 2: Perform NormCC Normalization
    run_normcc_normalization(metacc_folder)

    # Step 3: Run MetaCC Binning
    run_metacc_binning(metacc_folder, args.fasta, args.num_gene, args.seed)

    # Step 4: Run bin3C Binning
    run_bin3c_binning(bin3c_folder, args.fasta, args.seed)

    # Step 5: Run ImputeCC Binning
    # Ensure that contig_file and coverage file are provided
    if not args.contig_file:
        args.contig_file = os.path.join(metacc_folder, 'contig_info.csv')
    hic_matrix = os.path.join(metacc_folder, 'Normalized_contact_matrix.npz')
    run_imputecc_binning(imputecc_folder, args.contig_file, hic_matrix, args.fasta, args.num_gene)

    logger.info('All binning processes completed.')

if __name__ == "__main__":
    main()
