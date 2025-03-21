# hiczin_map.py

import numpy as np
import scipy.sparse as scisp
from math import log, exp
import logging

# package logger
logger = logging.getLogger(__name__)

class HiCzinMap:
    def __init__(self, path, contig_info, seq_map, norm_result, min_signal):
        '''
        perc: threshold of spurious contacts
        min_signal: minimum signal of acceptable contigs
        '''
        self.path = path
        self.seq_map = seq_map
        self.norm_result = norm_result
        self.signal = min_signal
        self.name = []
        self.site = []
        self.len = []
        self.cov = []
        self.tax = []

        for temp in contig_info:
            self.name.append(temp.name)
            self.site.append(temp.sites)
            self.len.append(temp.length)
            self.cov.append(temp.cov)
            self.tax.append(temp.tax)

        del contig_info

        self.name = np.array(self.name)
        self.site = np.array(self.site)
        self.len = np.array(self.len)
        self.cov = np.array(self.cov)
        self.tax = np.array(self.tax)
        
        # 避免覆盖度为零的情况
        self.cov[self.cov == 0] = np.min(self.cov[self.cov != 0])
        
        self.norm()
        del self.site, self.cov

    def norm(self):
        self.seq_map = self.seq_map.tocoo()
        _map_row = self.seq_map.row
        _map_col = self.seq_map.col
        _map_data = self.seq_map.data
        _map_coor = list(zip(_map_row, _map_col, _map_data))
        coeff = self.norm_result[0:4]  # 前四个回归系数

        self.seq_map = self.seq_map.tolil()
        self.seq_map = self.seq_map.astype(float)
        for x, y, d in _map_coor:

            s1 = self.site[x] if self.site[x] != 0 else 1
            s2 = self.site[y] if self.site[y] != 0 else 1
            s = (log(s1 * s2) - self.norm_result[5]) / self.norm_result[6]

            l1 = self.len[x]
            l2 = self.len[y]
            l = (log(l1 * l2) - self.norm_result[7]) / self.norm_result[8]

            c1 = self.cov[x]
            c2 = self.cov[y]
            c = (log(c1 * c2) - self.norm_result[9]) / self.norm_result[10]

            d_norm = d / exp(coeff[0] + coeff[1] * s + coeff[2] * l + coeff[3] * c)

            if d_norm > self.norm_result[4]:
                self.seq_map[x, y] = d_norm
            else:
                self.seq_map[x, y] = 0
        del _map_row, _map_col, _map_data, _map_coor