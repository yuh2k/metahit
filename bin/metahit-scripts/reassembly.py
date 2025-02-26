#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import sys
import os
import numpy as np
from scipy.stats import norm
import linecache
import random

script_directory = os.path.dirname(os.path.abspath(sys.argv[0]))

global_parser = argparse.ArgumentParser(add_help=False)
global_parser.add_argument('--bin', help='Binning result directory', required=True)
global_parser.add_argument('--hic1', help='Hi-C library forward fastq', required=True)
global_parser.add_argument('--hic2', help='Hi-C library reverse fastq', required=True)
global_parser.add_argument('--sg1', help='Shotgun forward fastq', required=True)
global_parser.add_argument('--sg2', help='Shotgun reverse fastq', required=True)
global_parser.add_argument('--bam', help='Hi-C bam that mapped to shotgun assembly', required=True)
global_parser.add_argument('--outdir', help='Output directory', required=True)
global_parser.add_argument('-p', '--metahit_path', help='Path to metahit folder', default=script_directory)
global_parser.add_argument('-t', help='threads', default='4')
global_parser.add_argument('-m', help='memory', default='24')

args = global_parser.parse_args()

os.mkdir(args.outdir)
fname  = os.path.join(args.outdir,"insert_size.txt")
print(fname)
insert_size_Cmd = f"samtools view {args.bam} | cut -f9 | awk '{{print sqrt($0^2)}}' > {fname}"
insert_size_execution = os.popen(insert_size_Cmd).read()
print(insert_size_execution)

def random_lines(filename):
    idxs = random.sample(range(10000), 10000)
    return [int(linecache.getline(filename, j)) for j in idxs if linecache.getline(filename, j).strip().isdigit()]

def EM (data):
    mu1_hat = 100
    sigma1_hat = 200
    pi1_hat=0.7
    mu2_hat = 800
    sigma2_hat = 500
    pi2_hat = 0.3

    difference = 10
    log_likelihood  = 10000
    while difference > 1e-1:
        gamma1 = pi1_hat * norm.pdf(data, mu1_hat, sigma1_hat)
        gamma2 = pi2_hat * norm.pdf(data, mu2_hat, sigma2_hat)
        total = gamma1 + gamma2
        if np.any(np.isnan(gamma1)):
            print("Array contains NaN values")
        gamma1 /= total
        gamma2 /= total

        mu1_hat = np.nansum(gamma1 * data) / np.sum(gamma1)
        mu2_hat = np.nansum(gamma2 * data) / np.sum(gamma2)
        sigma1_hat = np.sqrt(np.sum(gamma1 * (data - mu1_hat)**2) / np.sum(gamma1))
        sigma2_hat = np.sqrt(np.sum(gamma2 * (data - mu2_hat)**2) / np.sum(gamma2))
        pi1_hat = np.mean(gamma1)
        pi2_hat = np.mean(gamma2)

        difference = log_likelihood - np.sum(np.log(pi1_hat * norm.pdf(data, mu1_hat, sigma1_hat) + pi2_hat * norm.pdf(data, mu2_hat, sigma2_hat)))
        log_likelihood = np.sum(np.log(pi1_hat * norm.pdf(data, mu1_hat, sigma1_hat) + pi2_hat * norm.pdf(data, mu2_hat, sigma2_hat)))
    return mu1_hat,mu2_hat,sigma1_hat,sigma2_hat,pi1_hat,pi2_hat,difference,log_likelihood

with open(fname) as f:
    linecount = sum(1 for line in f)

data = []
for i in range(10):
    x = random_lines(fname)
    data = data + [z for z in x if z != 0]

data = np.array(data)
data = data[data < 1000]

mu1_hat,mu2_hat,sigma1_hat,sigma2_hat,pi1_hat,pi2_hat,difference,log_likelihood = EM(data)

print(difference, log_likelihood)
print("mu1_hat",mu1_hat)
print("mu2_hat",mu2_hat)
print("sigma1_hat",sigma1_hat)
print("sigma2_hat",sigma2_hat)
print("pi1_hat", pi1_hat)
print("pi2_hat", pi2_hat)

z_score = norm.ppf(0.05)
value = mu2_hat + z_score * sigma1_hat
print("value",value)

readname_sg_in_hic = os.path.join(args.outdir,"readname_sg_in_hic.txt")
sg_in_hic_forward = os.path.join(args.outdir,"sg_in_hic.forward.fastq")
sg_in_hic_reverse = os.path.join(args.outdir,"sg_in_hic.reverse.fastq")
new_hic_forward = os.path.join(args.outdir,"new_hic_forward.fastq")
new_hic_reverse = os.path.join(args.outdir,"new_hic_reverse.fastq")
new_sg_forward = os.path.join(args.outdir,"new_sg_forward.fastq")
new_sg_reverse = os.path.join(args.outdir,"new_sg_reverse.fastq")

filter_readname_cmd = f"samtools view {args.bam} | awk '($9>={value}*(-1) && $9<={value})' | cut -f1 | sort | uniq > {readname_sg_in_hic}"
filter_readname_forward_cmd = f"seqtk subseq {args.hic1} {readname_sg_in_hic} > {sg_in_hic_forward}"
filter_readname_reverse_cmd = f"seqtk subseq {args.hic2} {readname_sg_in_hic} > {sg_in_hic_reverse}"

hic_forward_cmd = f"cat {args.hic1} | paste - - - - | grep -v -F -f {readname_sg_in_hic} | tr '\\t' '\\n' > {new_hic_forward}"
hic_reverse_cmd = f"cat {args.hic2} | paste - - - - | grep -v -F -f {readname_sg_in_hic} | tr '\\t' '\\n' > {new_hic_reverse}"

cat_sg_forward_cmd = f"cat {sg_in_hic_forward} {args.sg1} > {new_sg_forward}"
cat_sg_reverse_cmd = f"cat {sg_in_hic_reverse} {args.sg2} > {new_sg_reverse}"

reassemble_cmd = f"bash {args.metahit_path}/bin/metahit-scripts/bin/modules/reassemble_bins.sh -b {args.bin} -o {args.outdir} -1 {args.sg1} -2 {args.sg2} -m {args.m} -t {args.t}"

output1 = os.popen(filter_readname_cmd).read()
print(output1)
output2 = os.popen(filter_readname_forward_cmd).read()
print(output2)
output3 = os.popen(filter_readname_reverse_cmd).read()
print(output3)
output4 = os.popen(hic_forward_cmd).read()
print(output4)
output5 = os.popen(hic_reverse_cmd).read()
print(output5)
output6 = os.popen(cat_sg_forward_cmd).read()
print(output6)
output7 = os.popen(cat_sg_reverse_cmd).read()
print(output7)
output8 = os.popen(reassemble_cmd).read()
print(output8)
