import scipy.sparse as sp

# Load the sparse matrix from the .npz file
matrix = sp.load_npz('../output/downstream/hic_contact_matrix.npz')


print(f"Matrix shape: {matrix.shape}")
print(f"Matrix type: {type(matrix)}")


dense_matrix = matrix.toarray()

# Display first 20 rows and cols
print(dense_matrix[:20, :20])  
