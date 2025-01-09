import argparse
import scipy.sparse as sp

# Setup argument parser
parser = argparse.ArgumentParser(description="Load and display a sparse matrix from an .npz file")
parser.add_argument('--file', required=True, help='Path to the .npz file')

args = parser.parse_args()

matrix = sp.load_npz(args.file)

print(f"Matrix shape: {matrix.shape}")
print(f"Matrix type: {type(matrix)}")

dense_matrix = matrix.toarray()
print(dense_matrix[:10, :10])

# python npz_view_tool.py --file output1/metacc/raw_contact_matrix_metacc.npz