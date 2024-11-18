#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
import logging

logger = logging.getLogger(__name__)

class NormCCMap:
    def __init__(self, path, contig_info, seq_map, norm_result, thres):
        self.path = path
        self.seq_map_raw = seq_map
        self.seq_map = seq_map
        self.norm_result = norm_result
        self.thres = thres
        
        # Check if contig_info is a DataFrame
        if not isinstance(contig_info, pd.DataFrame):
            raise ValueError("contig_info must be a pandas DataFrame")
        
        # Extract columns from the DataFrame
        self.name = contig_info['name'].tolist()
        self.site = contig_info['sites'].tolist()
        self.len = contig_info['length'].tolist()
        self.covcc = contig_info['covcc'].tolist()

        # Convert lists to numpy arrays for efficient processing
        self.name = np.array(self.name)
        self.site = np.array(self.site)
        self.len = np.array(self.len)
        self.covcc = np.array(self.covcc)

        # Perform normalization
        self.norm()

    def norm(self):
        self.seq_map = self.seq_map.tocoo()
        _map_row = self.seq_map.row
        _map_col = self.seq_map.col
        _map_data = self.seq_map.data
        _index = _map_row < _map_col
        _map_row = _map_row[_index]
        _map_col = _map_col[_index]
        _map_data = _map_data[_index]
        
        coeff = self.norm_result
        self.seq_map = self.seq_map.tolil()
        self.seq_map = self.seq_map.astype(float)

        mu_vector = [np.exp(coeff[0] + coeff[1] * np.log(site) + coeff[2] * np.log(length) + coeff[3] * np.log(covcc))
                     for site, length, covcc in zip(self.site, self.len, self.covcc)]
        scal = np.max(mu_vector)
        _norm_contact = []

        for i, j, d in zip(_map_row, _map_col, _map_data):
            d_norm = scal * d / np.sqrt(mu_vector[i] * mu_vector[j])
            _norm_contact.append(d_norm)
            self.seq_map[i, j] = d_norm
            self.seq_map[j, i] = d_norm

        logger.info('Eliminating systematic biases finished')

        # Remove spurious contacts based on threshold
        cutoffs = np.percentile(_norm_contact, self.thres * 100)
        for idx, val in enumerate(_norm_contact):
            if val < cutoffs:
                self.seq_map[_map_row[idx], _map_col[idx]] = 0
                self.seq_map[_map_col[idx], _map_row[idx]] = 0

        logger.info('Spurious contact detection finished')

        
        
        
        
class NormCCMap_LC:
    def __init__(self, path , contig_info , seq_map , norm_result , thres):
        '''
        perc: threshold of spurious contacts
        '''
        self.path = path
        self.seq_map_raw = seq_map
        self.seq_map = seq_map
        self.norm_result = norm_result
        self.thres = thres
        self.name = []
        self.len = []
        self.covcc = []

        for i in range(len(contig_info)):
            temp = contig_info[i]
            self.name.append(temp.name)
            self.len.append(temp.length)
            self.covcc.append(temp.covcc)

        del contig_info
        
        ####transfer the list to numpy array to do slicing#####
        self.name = np.array(self.name)
        self.len = np.array(self.len)
        self.covcc = np.array(self.covcc)
        
        #####Normalize raw contacts######
        self.norm()

        

    def norm(self):
        self.seq_map = self.seq_map.tocoo()
        _map_row = self.seq_map.row
        _map_col = self.seq_map.col
        _map_data = self.seq_map.data
        _index = _map_row<_map_col
        _map_row = _map_row[_index]
        _map_col = _map_col[_index]
        _map_data = _map_data[_index]
        
        _map_coor = list(zip(_map_row , _map_col , _map_data))
        coeff = self.norm_result
        
        self.seq_map = self.seq_map.tolil()
        self.seq_map = self.seq_map.astype(float)
        
        mu_vector = []
        for contig_feature in zip(self.len, self.covcc):
            mu_vector.append(exp(coeff[0] + coeff[1]*log(contig_feature[0])+ coeff[2]*log(contig_feature[1])))
        scal = np.max(mu_vector)
        _norm_contact = []
        
        for i in _map_coor:
            x = i[0]
            y = i[1]
            d = i[2]
            
            d_norm = scal*d/sqrt(mu_vector[x]*mu_vector[y])
            _norm_contact.append(d_norm)
            
            self.seq_map[x , y] = d_norm
            self.seq_map[y , x] = d_norm
            
        logger.info('Eliminating systematic biases finished')
        
        ########Remove spurious contacts###########
        cutoffs = np.percentile(_norm_contact , self.thres*100)
        count = 0
        for j in range(len(_norm_contact)):
            x = _map_row[j]
            y = _map_col[j]
            if _norm_contact[j] < cutoffs:
                self.seq_map[x , y] = 0
                self.seq_map[y , x] = 0
                count += 1
        logger.debug('{}% contacts have been removed with the cutoff {}'.format(round(100*count/len(_norm_contact)) , round(cutoffs,2)))
        logger.info('Spurious contact detection finished')
        
        del _map_row, _map_col, _map_data, _map_coor, _norm_contact, count



