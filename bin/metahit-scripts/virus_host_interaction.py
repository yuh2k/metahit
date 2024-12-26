#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 17:02:16 2024

@author: repair
"""
import argparse
import sys
import os
import numpy as np
from scipy.stats import norm


from MetaCC.Script.utils import load_object




script_directory = os.path.dirname(os.path.abspath(sys.argv[0]))

global_parser = argparse.ArgumentParser(add_help=False)
global_parser.add_argument('-V', '--version', default=False, action='store_true', help='Show the application version')
global_parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')

global_parser.add_argument('--BIN', help='binning result, a file that contain binning result')
global_parser.add_argument('--viral-contig', help='viral contig list')
global_parser.add_argument('--contact', help='contact matrix')
global_parser.add_argument('--OUTDIR', help='Output directory')
global_parser.add_argument('-p', help='metahit folder path')
global_parser.add_argument('-t', help='threads')
global_parser.add_argument('-m', help='memory')




args = global_parser.parse_args()

viral_contigs = []
host_contigs = []


with open(args.viral_contig) as vinputfile:
    for line in vinputfile:
        viral_contigs.append(line.strip())
        

with open(args.BIN) as hinputfile:
    for line in hinputfile:
        host_contigs.append(line.strip().split("\t"))
        
        



contact_matrix = load_object(args.contact)

seq_info = contact_matrix.seq_info

seq_map = contact_matrix.seq_map

result = []


for i in range(len(host_contigs)):
    for j in range(len(viral_contigs)):
        temp = []
        
        name_host = host_contigs[i][0]
        name_virus = viral_contigs[j]
        
        refid_host = seq_info[host_contigs[i][0]]["refid"]
        refid_virus = seq_info[viral_contigs[j]]["refid"]
                               
 
        temp.append(name_host)
        temp.append(name_virus)
        
        if refid_host < refid_virus: 
           temp.append(seq_map[refid_host,refid_virus])
           
        else:
            
           temp.append(seq_map[refid_virus,refid_host])
        result.append(temp)

print(result,file=args.OUTDIR)
        
        














