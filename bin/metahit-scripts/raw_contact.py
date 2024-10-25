#!/usr/bin/env python
# coding: utf-8

import os
import logging
from collections import OrderedDict, namedtuple
from collections.abc import Iterable
from Bio.Restriction import Restriction
import Bio.SeqIO as SeqIO
import numpy as np
import pysam
import scipy.sparse as scisp
import tqdm
from raw_utils import count_fasta_sequences, open_input
import argparse

# Configure logger to output debug information
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# Package logger
logger = logging.getLogger(__name__)

SeqInfo = namedtuple('SeqInfo', ['localid', 'refid', 'name', 'sites', 'length', 'covcc'])  # create a new class of tuple: SeqInfo

class SiteCounter(object):
    def __init__(self, enzyme_names, is_linear=True):
        """
        Simple class to count the total number of enzymatic cut sites for the given
        list of enzymes.
        :param enzyme_names: a list of enzyme names (proper case-sensitive spelling a la NEB)
        :param is_linear: Treat sequence as linear.
        """
        if isinstance(enzyme_names, str):
            enzyme_names = [enzyme_names]
        assert (isinstance(enzyme_names, Iterable) and
                not isinstance(enzyme_names, str)), 'enzyme_names must be a collection of names'
        self.enzymes = [getattr(Restriction, en) for en in enzyme_names]  # get the enzyme name from the standard module
        self.is_linear = is_linear

    def count_sites(self, seq):
        """
        Count the number of sites found in the given sequence, where sites from
        all specified enzymes are combined.
        :param seq: Bio.Seq object
        :return: the total number of sites
        """
        return sum(len(en.search(seq, self.is_linear)) for en in self.enzymes)


class Sparse2DAccumulator(object):
    # Create a 2D coo sparse matrix
    def __init__(self, N):
        self.shape = (N, N)
        self.mat = {}
        self.dtype = np.uint32

    def setitem(self, index, value):
        assert len(index) == 2 and index[0] >= 0 and index[1] >= 0, 'invalid index: {}'.format(index)
        self.mat[index] = value

    def getitem(self, index):
        if index in self.mat:
            return self.mat[index]
        else:
            return 0

    def get_coo(self, symm=True):
        """
        Create a COO format sparse representation of the accumulated values.
        :param symm: ensure matrix is symmetric on return
        :return: a scipy.coo_matrix sparse matrix
        """
        coords = [[], []]
        data = []
        m = self.mat
        for i, j in m.keys():
            coords[0].append(i)
            coords[1].append(j)
            data.append(m[i, j])

        m = scisp.coo_matrix((data, coords), shape=self.shape, dtype=self.dtype)

        if symm:
            m += scisp.tril(m.T, k=-1)

        return m.tocoo()


class ContactMatrix:
    def __init__(self, bam_file, enzymes, seq_file, path, min_mapq=30, min_len=1000, min_match=30, min_signal=2, coverage_file=None):
        # bam_file: alignment info of Hi-C library on contigs in bam
        # enzymes: name of restriction enzymes used in Hi-C experiments
        # seq_file: store the assembly contigs in fasta
        # path: output path
        # min_mapq: minimal mapping quality (default 30)
        # min_len: minimal length of contigs (default 1000bp)

        self.bam_file = bam_file
        self.enzymes = enzymes
        self.seq_file = seq_file
        self.path = path
        self.min_mapq_user = min_mapq
        self.min_len = min_len
        self.min_match_user = min_match
        self.min_signal = min_signal
        self.fasta_info = {}
        self.seq_info = []
        self.total_reads = None
        self.coverage_file = coverage_file

        # Parameters for metacc and bin3C
        self.min_mapq_metacc = 30
        self.min_match_metacc = 30
        self.min_mapq_bin3c = 60
        self.min_match_bin3c = 10

        # Ensure output directories exist
        os.makedirs(os.path.join(self.path, 'tmp'), exist_ok=True)

        logger.info('Reading fasta file...')
        with open_input(seq_file) as multi_fasta:
            # Prepare the site counter for the given experimental conditions
            site_counter = SiteCounter(enzymes, is_linear=True)
            # Get an estimate of sequences for progress
            fasta_count = count_fasta_sequences(seq_file)
            for seqrec in tqdm.tqdm(SeqIO.parse(multi_fasta, 'fasta'), total=fasta_count, desc='Analyzing contigs in reference fasta'):
                if len(seqrec) < self.min_len:
                    continue
                self.fasta_info[seqrec.id] = {'sites': site_counter.count_sites(seqrec.seq),
                                              'length': len(seqrec)}

        logger.debug('There are {} contigs in reference fasta'.format(fasta_count))

        # Now parse the header information from bam file
        with pysam.AlignmentFile(bam_file, 'rb') as bam:
            if 'SO' not in bam.header['HD'] or bam.header['HD']['SO'] != 'queryname':
                raise IOError('BAM file must be sorted by read name')

            logger.info('Filtering contigs by minimal length({})...'.format(self.min_len))
            ref_count = {'seq_missing': 0, 'too_short': 0}
            offset = 0
            localid = 0

            for n, (rname, rlen) in enumerate(zip(bam.references, bam.lengths)):
                if rlen < self.min_len:
                    ref_count['too_short'] += 1
                    continue

                try:
                    fa = self.fasta_info[rname]
                except KeyError:
                    logger.info('Contig "{}" was not present in reference fasta'.format(rname))
                    ref_count['seq_missing'] += 1
                    continue

                assert fa['length'] == rlen, 'Sequence lengths in {} do not agree: bam {} fasta {}'.format(rname, fa['length'], rlen)

                self.seq_info.append(SeqInfo(localid, n, rname, fa['sites'], rlen, 0))  # Initially set the covcc coverage to be zero
                localid += 1
                offset += rlen

            self.total_len = offset
            self.total_seq = localid

            del self.fasta_info

            if self.total_seq == 0:
                raise ImportError('No sequences in BAM file can be found in FASTA file')

            logger.debug('{} contigs missing and {} contigs are too short'.format(ref_count['seq_missing'], ref_count['too_short']))
            logger.debug('Accepted {} contigs covering {} bp'.format(self.total_seq, self.total_len))
            logger.info('Counting reads in bam file...')

            self.total_reads = bam.count(until_eof=True)
            logger.debug('BAM file contains {0} alignments'.format(self.total_reads))

            logger.info('Handling the alignments...')
            self._bin_map(bam)

        # Write contig_info.csv
        self.write_contig_info()
        self.save_contact_matrices()

    def _bin_map(self, bam):
        import tqdm

        def _matcher_user(r):
            if self.min_match_user > 0:
                return r.mapping_quality >= self.min_mapq_user and r.cigarstring is not None and r.cigartuples[0][0] == 0 and r.cigartuples[0][1] >= self.min_match_user
            else:
                return r.mapping_quality >= self.min_mapq_user

        def _matcher_metacc(r):
            if self.min_match_metacc > 0:
                return r.mapping_quality >= self.min_mapq_metacc and r.cigarstring is not None and r.cigartuples[0][0] == 0 and r.cigartuples[0][1] >= self.min_match_metacc
            else:
                return r.mapping_quality >= self.min_mapq_metacc

        def _matcher_bin3c(r):
            return True  # Accept all alignments

        def next_informative(_bam_iter, _pbar):
            while True:
                r = next(_bam_iter)
                _pbar.update()
                if not r.is_unmapped and not r.is_secondary and not r.is_supplementary:
                    return r

        _seq_map_user = Sparse2DAccumulator(self.total_seq)
        _seq_map_metacc = Sparse2DAccumulator(self.total_seq)
        _seq_map_bin3c = Sparse2DAccumulator(self.total_seq)

        with tqdm.tqdm(total=self.total_reads) as pbar:
            _idx = self.make_reverse_index('refid')
            bam.reset()
            bam_iter = bam.fetch(until_eof=True)

            counts = OrderedDict({
                'accepted map_different_contig pairs': 0,
                'accepted map_same_contig pairs': 0,
                'ref_excluded pairs': 0,
                'poor_match pairs': 0,
                'single read': 0
            })

            while True:
                try:
                    r1 = next_informative(bam_iter, pbar)
                    while True:
                        r2 = next_informative(bam_iter, pbar)
                        if r1.query_name == r2.query_name:
                            break
                        r1 = r2
                        counts['single read'] += 1
                except StopIteration:
                    break

                if r1.reference_id not in _idx or r2.reference_id not in _idx:
                    counts['ref_excluded pairs'] += 1
                    continue

                ix1 = _idx[r1.reference_id]
                ix2 = _idx[r2.reference_id]

                if ix2 < ix1:
                    ix1, ix2 = ix2, ix1

                ix = (ix1, ix2)

                # User-determined
                if _matcher_user(r1) and _matcher_user(r2):
                    if _seq_map_user.getitem(ix):
                        _seq_map_user.setitem(ix, _seq_map_user.getitem(ix) + 1)
                    else:
                        _seq_map_user.setitem(ix, 1)

                # MetaCC
                if _matcher_metacc(r1) and _matcher_metacc(r2):
                    if _seq_map_metacc.getitem(ix):
                        _seq_map_metacc.setitem(ix, _seq_map_metacc.getitem(ix) + 1)
                    else:
                        _seq_map_metacc.setitem(ix, 1)

                # Bin3C
                if _matcher_bin3c(r1) and _matcher_bin3c(r2):
                    if _seq_map_bin3c.getitem(ix):
                        _seq_map_bin3c.setitem(ix, _seq_map_bin3c.getitem(ix) + 1)
                    else:
                        _seq_map_bin3c.setitem(ix, 1)

            self.seq_map_user = _seq_map_user.get_coo()
            self.seq_map_metacc = _seq_map_metacc.get_coo()
            self.seq_map_bin3c = _seq_map_bin3c.get_coo()

            logger.debug('Pair accounting: {}'.format(counts))

    def make_reverse_index(self, field_name):
        rev_idx = {}
        for n, seq in enumerate(self.seq_info):
            fv = getattr(seq, field_name)
            if fv in rev_idx:
                raise RuntimeError('Field contains non-unique entries, a 1-1 mapping cannot be made')
            rev_idx[fv] = n
        return rev_idx

    def write_contig_info(self):
        contig_info_file = os.path.join(self.path, 'contig_info.csv')
        coverage_dict = {}

        if self.coverage_file and os.path.exists(self.coverage_file):
            logger.info(f'Coverage file found at {self.coverage_file}. Including coverage information in contig_info.csv.')
            with open(self.coverage_file, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) != 2:
                        continue
                    contig_name, coverage = parts
                    coverage_dict[contig_name] = float(coverage)
            has_coverage = True
        else:
            logger.info('Coverage file not provided or not found. Skipping coverage information in contig_info.csv.')
            has_coverage = False

        # Write contig_info.csv
        with open(contig_info_file, 'w') as f:
            for seq in self.seq_info:
                if has_coverage:
                    coverage = coverage_dict.get(seq.name, 0.0)
                    f.write(f'{seq.name},{seq.sites},{seq.length},{coverage}\n')
                else:
                    f.write(f'{seq.name},{seq.sites},{seq.length}\n')

        logger.info(f'contig_info.csv saved to {contig_info_file}')



    def save_contact_matrices(self):
        # Save the contact matrices as .npz files
        from scipy.sparse import save_npz

        matrices = {
            'user': self.seq_map_user,
            'metacc': self.seq_map_metacc,
            'bin3c': self.seq_map_bin3c
        }

        for key, matrix in matrices.items():
            matrix_file = os.path.join(self.path, f'contact_matrix_{key}.npz')
            save_npz(matrix_file, matrix)
            logger.info(f'Contact matrix "{key}" saved to {matrix_file}')

def main():
    parser = argparse.ArgumentParser(description="Generate raw contact matrix and contig info.")
    parser.add_argument('--bam', required=True, help='Path to the BAM file containing Hi-C reads')
    parser.add_argument('--fasta', required=True, help='Path to the FASTA file containing contig sequences')
    parser.add_argument('--out', required=True, help='Output directory to save the contact matrix and contig info')
    parser.add_argument('--enzymes', nargs='+', required=False, default=['HindIII'], help='List of enzymes')
    parser.add_argument('--min_mapq', type=int, default=30, help='Minimum MAPQ for user-determined method')
    parser.add_argument('--min_len', type=int, default=1000, help='Minimum contig length')
    parser.add_argument('--min_match', type=int, default=30, help='Minimum match length for user-determined method')
    parser.add_argument('--coverage', required=False, help='Path to the coverage.txt file')
    args = parser.parse_args()

    # Create ContactMatrix instance and process data
    cm = ContactMatrix(
        bam_file=args.bam,
        enzymes=args.enzymes,
        seq_file=args.fasta,
        path=args.out,
        min_mapq=args.min_mapq,
        min_len=args.min_len,
        min_match=args.min_match,
        coverage_file=args.coverage
    )


if __name__ == "__main__":
    main()
