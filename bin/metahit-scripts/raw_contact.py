#!/usr/bin/env python
# coding: utf-8

from collections import OrderedDict, namedtuple
from collections.abc import Iterable
from Bio.Restriction import Restriction
import Bio.SeqIO as SeqIO
import numpy as np
import pysam
import scipy.sparse as scisp
import tqdm
import os
from raw_utils import count_fasta_sequences, open_input
import logging
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
    def __init__(self, bam_file, enzymes, seq_file, path, min_mapq=30, min_len=1000, min_match=30, min_signal=2):
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
        self.min_mapq = min_mapq
        self.min_len = min_len
        self.min_match = min_match
        self.min_signal = min_signal
        self.fasta_info = {}
        self.seq_info = []
        self.seq_map = None
        self.total_reads = None

        # Ensure output directories exist
        os.makedirs(os.path.join(self.path, 'tmp'), exist_ok=True)

        logger.info('Reading fasta file...')
        with open_input(seq_file) as multi_fasta:
            # Prepare the site counter for the given experimental conditions
            site_counter = SiteCounter(enzymes, is_linear=True)
            # Get an estimate of sequences for progress
            fasta_count = count_fasta_sequences(seq_file)
            for seqrec in tqdm.tqdm(SeqIO.parse(multi_fasta, 'fasta'), total=fasta_count, desc='Analyzing contigs in reference fasta'):
                if len(seqrec) < min_len:
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
                if rlen < min_len:
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

        logger.info('Filtering contigs according to minimal signal({})...'.format(self.min_signal))
        contig_id = self.max_offdiag()
        logger.debug('{} contigs remain'.format(len(contig_id)))

        self.seq_map = self.seq_map.tolil()
        seq_temp = []  # Temporarily store the sequence information
        for i, idn in enumerate(contig_id):
            seq = self.seq_info[idn]
            assert seq.localid == idn, 'The local index does not match the contact matrix index'
            seq_temp.append(SeqInfo(i, seq.refid, seq.name, seq.sites, seq.length, self.seq_map[idn, idn]))

        self.seq_info = seq_temp
        del seq_temp

        self.seq_map = self.seq_map.tocsr()
        self.seq_map = self.seq_map[contig_id, :]
        self.seq_map = self.seq_map.tocsc()
        self.seq_map = self.seq_map[:, contig_id]
        self.seq_map = self.seq_map.tocoo()
        del contig_id

        assert self.seq_map.shape[0] == len(self.seq_info), 'Filter error'



    def _bin_map(self, bam):
        import tqdm

        def _simple_match(r):
            return r.mapping_quality >= self.min_mapq

        def _strong_match(r):
            if r.mapping_quality < self.min_mapq or r.cigarstring is None:
                return False
            cig = r.cigartuples[-1] if r.is_reverse else r.cigartuples[0]
            return cig[0] == 0 and cig[1] >= self.min_match

        _matcher = _strong_match if self.min_match else _simple_match

        def next_informative(_bam_iter, _pbar):
            while True:
                r = next(_bam_iter)
                _pbar.update()
                if not r.is_unmapped and not r.is_secondary and not r.is_supplementary:
                    return r

        _seq_map = Sparse2DAccumulator(self.total_seq)

        with tqdm.tqdm(total=self.total_reads) as pbar:
            _mapq = self.min_mapq
            _idx = self.make_reverse_index('refid')
            bam.reset()
            bam_iter = bam.fetch(until_eof=True)
            self.index1 = 0

            counts = OrderedDict({
                'accepted map_different_contig pairs': 0,
                'accepted map_same_contig pairs': 0,
                'ref_excluded pairs': 0,
                'poor_match pairs': 0,
                'single read': 0
            })

            while True:
                self.index1 += 1
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

                if not _matcher(r1) or not _matcher(r2):
                    counts['poor_match pairs'] += 1
                    continue

                if r1.reference_id == r2.reference_id:
                    counts['accepted map_same_contig pairs'] += 1
                else:
                    counts['accepted map_different_contig pairs'] += 1

                ix1 = _idx[r1.reference_id]
                ix2 = _idx[r2.reference_id]

                if ix2 < ix1:
                    ix1, ix2 = ix2, ix1

                ix = (ix1, ix2)
                if _seq_map.getitem(ix):
                    _seq_map.setitem(ix, _seq_map.getitem(ix) + 1)
                else:
                    _seq_map.setitem(ix, 1)

        self.seq_map = _seq_map.get_coo()
        del _seq_map, r1, r2, _idx

        logger.debug('Pair accounting: {}'.format(counts))

    def make_reverse_index(self, field_name):
        rev_idx = {}
        for n, seq in enumerate(self.seq_info):
            fv = getattr(seq, field_name)
            if fv in rev_idx:
                raise RuntimeError('Field contains non-unique entries, a 1-1 mapping cannot be made')
            rev_idx[fv] = n
        return rev_idx



    def max_offdiag(self):
        _m = self.seq_map
        assert scisp.isspmatrix(_m), 'Input matrix is not a scipy.sparse object'
        _m = _m.tolil(True)
        _diag = _m.tocsr().diagonal()
        _m.setdiag(0)
        _sig = np.asarray(_m.tocsr().max(axis=0).todense()).ravel()
        _contig_id = [i for i in range(_m.shape[0]) if _sig[i] >= self.min_signal and _diag[i] > 0 and self.seq_info[i].sites > 0]
        del _m
        return _contig_id


# Configure logger to output debug information
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# package logger
logger = logging.getLogger(__name__)

# offset is the enumeration of the length
# localid is the index of the list
# refid is the index of the fasta file, which is a global index
SeqInfo = namedtuple('SeqInfo', ['localid', 'refid', 'name', 'sites', 'length', 'covcc'])  # create a new class of tuple: SeqInfo


