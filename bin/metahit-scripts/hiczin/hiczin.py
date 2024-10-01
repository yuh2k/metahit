#########The structure of the main script is modified from bin3C########
from raw_contact import ContactMatrix
from hiczin_contact import HiCzinMap
from exceptions import ApplicationException
from hiczin_utils import load_object, save_object, make_dir
import scipy.sparse as scisp
import logging
import sys
import argparse
import os
import time
import warnings

##Ignore the warning information of package deprecation##
warnings.filterwarnings("ignore")

__version__ = '0.1.0, released at 12/2020'

if __name__ == '__main__':
    
    def mk_version():
        return 'HiCzin v{}'.format(__version__)

    def out_name(base, suffix):
        return '{}{}'.format(base, suffix)

    def ifelse(arg, default):
        if arg is None:
            return default
        else:
            return arg

    runtime_defaults = {
        'min_len': 1000,
        'min_signal': 2,
        'min_mapq': 30,
        'min_match': 30,
        'thres': 0.05
    }

    global_parser = argparse.ArgumentParser(add_help=False)
    global_parser.add_argument('-V', '--version', default=False, action='store_true', help='Show the application version')
    global_parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    global_parser.add_argument('--cover', default=False, action='store_true', help='Cover existing files')
    global_parser.add_argument('--log', help='Log file path [OUTDIR/hiczin.log]')


    parser = argparse.ArgumentParser(description='HiCzin: Normalizing metagenomic Hi-C data and detecting spurious contacts using zero-inflated negative binomial regression')

    subparsers = parser.add_subparsers(title='commands', dest='command', description='Valid commands',
                                       help='choose an analysis stage for further options')

    cmd_pl = subparsers.add_parser('norm', parents=[global_parser],
                                      description='Normalize contacts and do the binning.')

    cmd_test = subparsers.add_parser('test', parents=[global_parser],
                                        description='pipeline testing.')

    '''
    pipeline subparser input
    '''
    cmd_pl.add_argument('--min-len', type=int,
                               help='Minimum acceptable reference length [1000]')
    cmd_pl.add_argument('--min-signal', type=int, help='Minimum acceptable signal [2]')
    cmd_pl.add_argument('--min-mapq', type=int,
                           help='Minimum acceptable mapping quality [30]')
    cmd_pl.add_argument('--min-match', type=int,
                           help='Accepted alignments must being N matches [30]')
    cmd_pl.add_argument('-e', '--enzyme', metavar='NEB_NAME', action='append',
                           help='Case-sensitive enzyme name. Use multiple times for multiple enzymes')
    cmd_pl.add_argument('--thres', type=float,
                           help='acceptable fraction of incorrectly identified valid contacts [0.05]')
    cmd_pl.add_argument('FASTA', help='Reference fasta sequence')
    cmd_pl.add_argument('BAM', help='Input bam file in query order')
    cmd_pl.add_argument('TAX', help='Contig labels from TAXAssign')
    cmd_pl.add_argument('COV', help='Coverage information of assembled contigs')
    cmd_pl.add_argument('OUTDIR', help='Output directory')
         

    '''
    Testing of HiCzin software
    '''
    cmd_test.add_argument('OUTDIR', help='Output directory of testing results')  

    args = parser.parse_args()

    if args.version:
        print(mk_version())
        sys.exit(0)

    try:
        make_dir(args.OUTDIR, args.cover)
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
    if args.verbose:
        ch.setLevel(logging.DEBUG)
    else:
        ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    root.addHandler(ch)

    # File log listens to all levels from root
    if args.log is not None:
        log_path = args.log
    else:
        log_path = os.path.join(args.OUTDIR, 'hiczin.log')
    fh = logging.FileHandler(log_path, mode='a')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    root.addHandler(fh)

    # Add some environmental details
    logger.debug(mk_version())
    logger.debug(sys.version.replace('\n', ' '))
    logger.debug('Command line: {}'.format(' '.join(sys.argv)))

    try:
        start_time = time.time()
        if args.command == 'norm':
            if args.enzyme is not None: 
            # Create a contact map for analysis by HiCzin
                cm = ContactMatrix(args.BAM,
                                args.enzyme,
                                args.FASTA,
                                args.TAX,
                                args.COV,
                                args.OUTDIR,
                                min_mapq=ifelse(args.min_mapq, runtime_defaults['min_mapq']),
                                min_len=ifelse(args.min_len, runtime_defaults['min_len']),
                                min_match=ifelse(args.min_match, runtime_defaults['min_match']),
                                min_signal=ifelse(args.min_signal, runtime_defaults['min_signal']))
                
                if cm.is_empty():
                    logger.info('Stopping as the contact matrix is empty')
                    sys.exit(1)

                logger.info('Raw contact matrix construction finished')

                logger.info('Normalizing raw contacts by HiCzin...')
                
                from rpy2 import robjects
                r = robjects.r
                r.source('HiCzin.R')
                contig_file = os.path.join(args.OUTDIR ,'contig_info.csv')
                valid_file = os.path.join(args.OUTDIR ,'valid_contact.csv')
                thres = ifelse(args.thres, runtime_defaults['thres'])
                norm_result = r.HiCzin(contig_file , valid_file , thres)
                
                hzmap = HiCzinMap(args.OUTDIR,
                                cm.seq_info,
                                cm.seq_map,
                                norm_result,
                                ifelse(args.min_signal, runtime_defaults['min_signal']))
                               

                

            else:
                cm = ContactMatrix_LC(args.BAM,
                                args.FASTA,
                                args.TAX,
                                args.COV,
                                args.OUTDIR,
                                min_mapq=ifelse(args.min_mapq, runtime_defaults['min_mapq']),
                                min_len=ifelse(args.min_len, runtime_defaults['min_len']),
                                min_match=ifelse(args.min_match, runtime_defaults['min_match']),
                                min_signal=ifelse(args.min_signal, runtime_defaults['min_signal']))

                if cm.is_empty():
                    logger.info('Stopping as the contact matrix is empty')
                    sys.exit(1)

                logger.info('Raw contact matrix construction finished')

                logger.info('Normalizing by HiCzin_LC because no enzymes are input...')
            
                from rpy2 import robjects
                r = robjects.r
                r.source('HiCzin.R')
                contig_file = os.path.join(args.OUTDIR ,'contig_info.csv')
                valid_file = os.path.join(args.OUTDIR ,'valid_contact.csv')
                thres = ifelse(args.thres, runtime_defaults['thres'])
                norm_result = r.HiCzin(contig_file , valid_file , thres)
            
                hzmap = HiCzinMap_LC(args.OUTDIR,
                                cm.seq_info,
                                cm.seq_map,
                                norm_result,
                                ifelse(args.min_signal, runtime_defaults['min_signal']))

            logger.info('Normalization finished')
            

            end_time = time.time()
            logger.info('HiCzin consumes {} seconds in total'.format(str(end_time-start_time)))

            logger.info('Saving Normalized Hi-C contact maps')
            scisp.save_npz(os.path.join(args.OUTDIR, 'Normalized_contact_matrix.npz'), hzmap.seq_map.tocsr())
            save_object(os.path.join(args.OUTDIR, 'HiCzin_normalized_contact'), hzmap)



        if args.command == 'test':
            logger.info('Begin to test the HiCzin software...')
            ENZ = 'HindIII'
            FASTA = 'test/contig_assembly_test.fa'
            BAM = 'test/align_test.bam'
            TAX = 'test/contig_tax_test.csv'
            COV = 'test/coverage_test.txt'
            logger.info('Begin to test the contact map construction section...')
            cm = ContactMatrix(BAM,
                                ENZ,
                                FASTA,
                                TAX,
                                COV,
                                args.OUTDIR,
                                min_mapq=runtime_defaults['min_mapq'],
                                min_len=runtime_defaults['min_len'],
                                min_match=runtime_defaults['min_match'],
                                min_signal=0)

            logger.info('Contact map construction section works!')
            
            logger.info('Begin to test the normalization section...')
            logger.info('Normalizing by HiCzin...')

            
            from rpy2 import robjects
            r = robjects.r
            r.source('HiCzin.R')
            contig_file = 'test/contig_info_test.csv'
            valid_file = 'test/must_contact_test.csv'
            thres = runtime_defaults['thres']
            norm_result = r.HiCzin(contig_file , valid_file , thres)
            
            hzmap = HiCzinMap(args.OUTDIR,
                            cm.seq_info,
                            cm.seq_map,
                            norm_result,
                            0)
            scisp.save_npz(os.path.join(args.OUTDIR, 'Normalized_contact_matrix.npz'), hzmap.seq_map.tocsr())
            save_object(os.path.join(args.OUTDIR, 'HiCzin_normalized_contact'), hzmap)
            logger.info('Normalization section works!')
            logger.info('Testing finished!')




    except ApplicationException:
        logger.error('ApplicationException Error')
        sys.exit(1)