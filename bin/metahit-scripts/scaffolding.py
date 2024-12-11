#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 21:38:59 2024

@author: repair
"""

#######Import python packages
import argparse
import logging
import sys
import os

class ApplicationException(Exception):
    def __init__(self, message):
        super(ApplicationException, self).__init__(message)


if __name__ == '__main__':
    
    
    def make_dir(path, exist_ok=True):
        """
        Convenience method for making directories with a standard logic.
        An exception is raised when the specified path exists and is not a directory.
        :param path: target path to create
        :param exist_ok: if true, an existing directory is ok. Existing files will still cause an exception
        """
        if not os.path.exists(path):
            os.mkdir(path)
        elif not exist_ok:
            raise IOError('output directory already exists!')
        elif os.path.isfile(path):
            raise IOError('output path already exists and is a file!')
            
    
    global_parser = argparse.ArgumentParser(add_help=False)
    global_parser.add_argument('--fasta', help='Reference fasta sequence of bin for scaffolding')
    global_parser.add_argument('--bam', help='Input bam file of Hi-C reads mapped to shotgun assembly')
    global_parser.add_argument('--enzyme', help='list of restriction enzymes, separated by comma')
    global_parser.add_argument('--outdir', help='Output directory path')
    global_parser.add_argument('--hic1', help='Hi-C library forward fastq')
    global_parser.add_argument('--hic2', help='Hi-C library reverse fastq')
    global_parser.add_argument('-t', help='threads')
    global_parser.add_argument('-m', help='memory')
    global_parser.add_argument('-r', help='resolution of resulting graph, default is 10kb')
    args = global_parser.parse_args()
    
    try:
        make_dir(args.outdir)
    except IOError:
        print('Error: cannot find out directory or the directory already exists')
        sys.exit(1)
    

    
   
    
    logging.captureWarnings(True)
    logger = logging.getLogger('main')

    # root log listens to everything
    root = logging.getLogger('')
    root.setLevel(logging.DEBUG)

    # log message format
    formatter = logging.Formatter(fmt='%(levelname)-8s | %(asctime)s | %(name)7s | %(message)s')

    # Runtime console listens to INFO by default
    ch = logging.StreamHandler()
    
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    root.addHandler(ch)

    # File log listens to all levels from root
    log_path = os.path.join(args.outdir, 'scaffolding.log')
    fh = logging.FileHandler(log_path, mode='a')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    root.addHandler(fh)

    # Add some environmental details
    logger.debug(sys.version.replace('\n', ' '))
    logger.debug('Command line: {}'.format(' '.join(sys.argv)))
    

    

    try:
        
        yahs_path=os.path.join(args.outdir, "yahs")
 
        hicstuff_pipe_path=os.path.join(args.outdir, "hicstuff_pipeline")

        hicstuff_view_path=os.path.join(args.outdir, "hicstuff_view")
        
        
        yahs_cmd = "yahs " + args.fasta + " " + args.bam + " -o "+ yahs_path+"/final"
        logger.info('Running command ' + yahs_cmd)
        output = os.popen(yahs_cmd).read()
        logger.info('OUTPUT IS ' + output)
        
        hicstuff_pipe_cmd = f"./bin/metahit-scripts/hicstuff/hicstuff.py pipeline -a bwa -Df -g {yahs_path}/final_scaffolds_final.fa -e {args.enzyme} -o {hicstuff_pipe_path} -t {args.t} {args.hic1} {args.hic2}"
        
        hicstuff_view_cmd = f"./bin/metahit-scripts/hicstuff/hicstuff.py view -f {hicstuff_pipe_path}/fragments_list.txt -n -b  -o {hicstuff_view_path}/scaffolding.png {hicstuff_pipe_path}/abs_fragments_contacts_weighted.cool"
        
        output1 = os.popen(hicstuff_pipe_cmd).read()
        logger.info(output1)
        
        output2 = os.popen(hicstuff_view_cmd).read()
        logger.info(output2)
        
        
    except ApplicationException:
        logger.error('ApplicationException Error')
        sys.exit(1)
    
    
    
    
    
        
    
    