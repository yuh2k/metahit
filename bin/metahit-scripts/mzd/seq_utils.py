import bz2
import pickle
import gzip
import json
import io
import yaml
from Bio.Restriction import Restriction
from difflib import SequenceMatcher
from collections.abc import Iterable, Mapping
from scipy.stats import mstats
import Bio.SeqIO as SeqIO
import heapq
import numpy as np
import networkx as nx
import os
import subprocess
import tempfile
import uuid
import logging

logger = logging.getLogger(__name__)


def count_fasta_sequences(file_name):
    """
    Estimate the number of fasta sequences in a file by counting headers. Decompression is automatically attempted
    for files ending in .gz. Counting and decompression is by way of subprocess calls to grep and gzip. Uncompressed
    files are also handled. This is faster than parsing a file with BioPython.
    """
    if file_name.endswith('.gz'):
        proc_uncomp = subprocess.Popen(['gzip', '-cd', file_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc_read = subprocess.Popen(['grep', r'^>'], stdin=proc_uncomp.stdout, stdout=subprocess.PIPE)
    else:
        proc_read = subprocess.Popen(['grep', r'^>', file_name], stdout=subprocess.PIPE)

    n = 0
    for _ in proc_read.stdout:
        n += 1
    return n


class IndexedFasta(Mapping):
    """
    Provides indexed access to a sequence file in Fasta format which can be compressed using bgzip.
    """

    def __init__(self, fasta_file, tmp_path=None):
        if tmp_path is None:
            tmp_path = tempfile.gettempdir()
        elif not os.path.exists(tmp_path):
            raise IOError(f'specified temporary path [{tmp_path}] does not exist')
        self._tmp_file = os.path.join(tmp_path, f'{uuid.uuid4()}.seqdb')
        self._fasta_file = fasta_file
        self._index = SeqIO.index_db(self._tmp_file, self._fasta_file, 'fasta')

    def __getitem__(self, _id):
        return self._index[_id]

    def __iter__(self):
        return iter(self._index)

    def __len__(self):
        return len(self._index)

    def close(self):
        if self._index:
            self._index.close()
        if os.path.exists(self._tmp_file):
            os.unlink(self._tmp_file)


class SiteCounter:

    def __init__(self, enzyme_names, tip_size=None, is_linear=True):
        """
        Simple class to count the total number of enzymatic cut sites for the given list of enzymes.
        """
        if isinstance(enzyme_names, str):
            enzyme_names = [enzyme_names]
        assert isinstance(enzyme_names, Iterable) and not isinstance(enzyme_names, str), 'enzyme_names must be a collection of names'
        self.enzymes = [SiteCounter._get_enzyme_instance(en) for en in enzyme_names]
        self.is_linear = is_linear
        self.tip_size = tip_size

    @staticmethod
    def _get_enzyme_instance(enz_name):
        try:
            return getattr(Restriction, enz_name)
        except AttributeError:
            similar = [a for a in dir(Restriction) if SequenceMatcher(None, enz_name.lower(), a.lower()).ratio() >= 0.8]
            raise UnknownEnzymeException(enz_name, similar)

    def _count(self, seq):
        return sum(len(en.search(seq, self.is_linear)) for en in self.enzymes)

    def count_sites(self, seq):
        """
        Count the number of sites found in the given sequence, where sites from all specified enzymes are combined.
        """
        if self.tip_size:
            seq_len = len(seq)
            if seq_len < 2 * self.tip_size:
                l_tip = seq[:seq_len // 2]
                r_tip = seq[-seq_len // 2:]
            else:
                l_tip = seq[:self.tip_size]
                r_tip = seq[-self.tip_size:]
            sites = [self._count(l_tip), self._count(r_tip)]
        else:
            sites = self._count(seq)

        return sites


class SequenceAnalyzer:

    COV_TYPE = np.dtype([('index', int), ('status', np.bool_), ('node', float),
                         ('local', float), ('fold', float)])

    @staticmethod
    def read_report(file_name):
        return yaml.safe_load(open(file_name, 'r'))

    def __init__(self, seq_map, seq_report, seq_info, tip_size):
        self.seq_map = seq_map
        self.seq_report = seq_report
        self.seq_info = seq_info
        self.tip_size = tip_size

    def _contact_graph(self):
        g = nx.Graph()
        n_seq = len(self.seq_info)

        for i in range(n_seq):
            si = self.seq_info[i]
            d = self.seq_report['seq_info'][si.name]
            if self.tip_size:
                g.add_node(i, _id=si.name, _cov=float(d['coverage']), _sites=d['sites'], _len=int(d['length']))
            else:
                g.add_node(i, _id=si.name, _cov=float(d['coverage']), _sites=int(d['sites']), _len=int(d['length']))

        _m = self.seq_map
        if self.tip_size:
            _m = _m.sum(axis=(2, 3))

        for i in range(n_seq):
            for j in range(i, n_seq):
                if _m[i, j] > 0:
                    g.add_edge(i, j, weight=float(_m[i, j]))

        return g

    @staticmethod
    def _nlargest(g, u, n, k=0, local_set=None):
        """
        Build a list of nodes of length n within a radius k hops of node u in graph g, which have the largest weight.
        """
        if not local_set:
            local_set = set()

        neighbors = [v[0] for v in heapq.nlargest(n + 1, g[u].items(), key=lambda x: x[1]['weight'])]
        local_set.update(neighbors)
        if k > 0:
            for v in neighbors:
                if v == u:
                    continue
                SequenceAnalyzer._nlargest(g, v, n, k - 1, local_set)

        return list(local_set)

    def report_degenerates(self, fold_max, min_len=0):
        """
        Report nodes in the graph whose coverage exceeds a threshold.
        """
        g = self._contact_graph()
        degens = []

        for u in g.nodes():
            if g.nodes[u]['_len'] < min_len or g.degree(u) == 0:
                continue

            signif_local_nodes = SequenceAnalyzer._nlargest(g, u, 4, 1)
            local_mean_cov = mstats.gmean(np.array([g.nodes[v]['_cov'] for v in signif_local_nodes]))
            fold_vs_local = g.nodes[u]['_cov'] / float(local_mean_cov)

            is_degen = fold_vs_local > fold_max

            degens.append((u, is_degen, g.nodes[u]['_cov'], local_mean_cov, fold_vs_local))

        degens = np.array(degens, dtype=SequenceAnalyzer.COV_TYPE)

        if len(degens) == 0:
            logger.debug('No degenerate sequences found')
        else:
            logger.debug('Degenerate sequence report')
            for di in degens[degens['status']]:
                logger.debug(di)

        return degens
