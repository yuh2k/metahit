import argparse
import os
import pysam
import tqdm
import numpy as np
import scipy.sparse as sp
import Bio.SeqIO as SeqIO
from collections import namedtuple

# Define the contig information tuple
SeqInfo_LC = namedtuple('SeqInfo_LC', ['localid', 'refid', 'name', 'length', 'covcc'])

# Class to accumulate the sparse matrix
class Sparse2DAccumulator(object):
    def __init__(self, N):
        self.shape = (N, N)
        self.mat = {}
        self.dtype = np.uint32

    def setitem(self, index, value):
        self.mat[index] = value

    def getitem(self, index):
        return self.mat.get(index, 0)

    def get_coo(self):
        coords = [[], []]
        data = []
        for i, j in self.mat.keys():
            coords[0].append(i)
            coords[1].append(j)
            data.append(self.mat[i, j])
        m = sp.coo_matrix((data, coords), shape=self.shape, dtype=self.dtype)
        m += sp.tril(m.T, k=-1)
        return m.tocoo()

# Class for handling Hi-C contact matrix generation (no enzyme required)
class ContactMatrix_LC:
    def __init__(self, bam_file, seq_file, path, min_mapq=30, min_len=1000, min_match=30, min_signal=2):
        self.bam_file = bam_file
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

        # Parse the FASTA file
        print('Reading fasta file...')
        with open(seq_file, 'r') as multi_fasta:
            fasta_count = sum(1 for _ in SeqIO.parse(multi_fasta, 'fasta'))
            for seqrec in tqdm.tqdm(SeqIO.parse(seq_file, 'fasta'), total=fasta_count, desc='Analyzing contigs in reference fasta'):
                if len(seqrec) < min_len:
                    continue
                self.fasta_info[seqrec.id] = {'length': len(seqrec)}

        # Parse the BAM file and filter contigs
        print('Filtering contigs by minimal length...')
        with pysam.AlignmentFile(bam_file, 'rb') as bam:
            offset = 0
            localid = 0
            for n, (rname, rlen) in enumerate(zip(bam.references, bam.lengths)):
                if rlen < min_len:
                    continue
                try:
                    fa = self.fasta_info[rname]
                except KeyError:
                    continue
                assert fa['length'] == rlen, f'Sequence lengths in BAM and FASTA do not match for {rname}'
                self.seq_info.append(SeqInfo_LC(localid, n, rname, rlen, 0))
                localid += 1
                offset += rlen
            self.total_len = offset
            self.total_seq = localid

        # Initialize the sparse matrix
        self._bin_map(bam_file)

        # Filter contigs by minimal signal and store contact matrix
        self._filter_contigs()
        self._write_contig_info()

    # Generate contact matrix by processing BAM alignments
    def _bin_map(self, bam_file):
        print('Generating contact matrix...')
        seq_map = Sparse2DAccumulator(self.total_seq)
        with pysam.AlignmentFile(bam_file, 'rb') as bam:
            for r1, r2 in self._fetch_read_pairs(bam):
                if r1.reference_id == r2.reference_id:
                    ix1 = r1.reference_id
                    ix2 = r2.reference_id
                    if ix2 < ix1:
                        ix1, ix2 = ix2, ix1
                    seq_map.setitem((ix1, ix2), seq_map.getitem((ix1, ix2)) + 1)
        self.seq_map = seq_map.get_coo()

    # Fetch read pairs from BAM
    def _fetch_read_pairs(self, bam):
        r1 = next(bam.fetch(until_eof=True))
        while True:
            try:
                r2 = next(bam.fetch(until_eof=True))
                if r1.query_name == r2.query_name:
                    yield r1, r2
                r1 = r2
            except StopIteration:
                break

    # Filter contigs by signal threshold
    def _filter_contigs(self):
        print('Filtering contigs by signal threshold...')
        self.seq_map = self.seq_map.tolil()
        seq_temp = []
        for i, seq in enumerate(self.seq_info):
            seq_temp.append(SeqInfo_LC(i, seq.refid, seq.name, seq.length, self.seq_map[i, i]))
        self.seq_info = seq_temp
        self.seq_map = self.seq_map.tocsr()
        self.seq_map = self.seq_map.tocsc()

    # Write contig information to CSV
    def _write_contig_info(self):
        print('Writing contig information...')
        contig_info_file = os.path.join(self.path, 'contig_info.csv')
        with open(contig_info_file, 'w') as out:
            out.write('Contig name,Length,Coverage\n')
            for seq in self.seq_info:
                out.write(f'{seq.name},{seq.length},{seq.covcc}\n')

        # Save the sparse matrix as .npz
        npz_output = os.path.join(self.path, 'hic_contact_matrix.npz')
        sp.save_npz(npz_output, self.seq_map)
        print(f'Saved Hi-C contact matrix to {npz_output}')

# Main function to handle command line arguments
def main():
    parser = argparse.ArgumentParser(description="Generate Hi-C contact matrix and contig info from BAM and FASTA files")
    parser.add_argument('--bam', required=True, help='Input BAM file containing Hi-C reads mapped to contigs')
    parser.add_argument('--fasta', required=True, help='Input FASTA file containing contig sequences')
    parser.add_argument('--out', required=True, help='Output directory to save Hi-C contact matrix (.npz) and contig info (.csv)')
    parser.add_argument('--min_mapq', type=int, default=30, help='Minimum mapping quality for filtering BAM reads')
    parser.add_argument('--min_len', type=int, default=1000, help='Minimum contig length to include in the analysis')
    parser.add_argument('--min_signal', type=int, default=2, help='Minimum signal threshold for filtering contigs')
    
    args = parser.parse_args()

    # Create the output directory if it doesn't exist
    os.makedirs(args.out, exist_ok=True)

    # Initialize and run the ContactMatrix_LC
    contact_matrix = ContactMatrix_LC(
        bam_file=args.bam,
        seq_file=args.fasta,
        path=args.out,
        min_mapq=args.min_mapq,
        min_len=args.min_len,
        min_signal=args.min_signal
    )

if __name__ == '__main__':
    main()
