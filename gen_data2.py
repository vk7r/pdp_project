import os
import numpy as np
from scipy.sparse import coo_matrix

def generate_sparse_coo_matrix(n, density):
    nnz = int(n * n * density)
    nnz = max(nnz, n)  # ensure at least some data
    
    # Random row and col indices for nonzeros
    row = np.random.randint(0, n, size=nnz)
    col = np.random.randint(0, n, size=nnz)
    data = np.random.uniform(1, 10, size=nnz)
    
    matrix = coo_matrix((data, (row, col)), shape=(n, n))
    matrix.sum_duplicates()  # combine duplicates
    
    return matrix

def save_coo_matrix_market(matrix, filename):
    with open(filename, 'w') as f:
        f.write('%%MatrixMarket matrix coordinate real general\n')
        f.write(f"{matrix.shape[0]} {matrix.shape[1]} {matrix.nnz}\n")
        for r, c, v in zip(matrix.row, matrix.col, matrix.data):
            f.write(f"{r+1} {c+1} {v:.6g}\n")

def main():
    os.makedirs('coo_input_data', exist_ok=True)
    
    # Sizes and densities (density smaller for larger matrices)
    configs = [
        (10, 0.3),       # small, denser
        (100, 0.1),      # medium size, moderate density
        (10_000, 0.001), # large, very sparse
        (100_000, 0.00001)  # huge, extremely sparse
    ]
    
    for n, density in configs:
        print(f"Generating matrix {n}x{n} with density {density}")
        matrix = generate_sparse_coo_matrix(n, density)
        filename = f"coo_input_data/matrix_{n}x{n}_d{int(density*100)}.mtx"
        save_coo_matrix_market(matrix, filename)
        print(f"Saved {filename} with {matrix.nnz} nonzeros")

if __name__ == "__main__":
    main()
