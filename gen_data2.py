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

    # Define sizes
    sizes = [100, 1000, 10_000, 100_000]

    # Define three levels of density: sparse, medium, dense
    density_levels = [0.001, 0.01, 0.1]  # Adjust as needed for your case

    for n in sizes:
        for density in density_levels:
            print(f"Generating matrix {n}x{n} with density {density}")
            matrix = generate_sparse_coo_matrix(n, density)
            density_str = f"{density:.5f}".replace('.', '_')
            filename = f"new_coo_input_data/matrix_{n}x{n}_d{density_str}.mtx"
            save_coo_matrix_market(matrix, filename)
            print(f"Saved {filename} with {matrix.nnz} nonzeros")

if __name__ == "__main__":
    main()

# def main():
#     os.makedirs('coo_input_data', exist_ok=True)
    
#     # Sizes and densities (density smaller for larger matrices)
#     configs = [
#         (100, 0.3),
#         (1000, 0.1),
#         (10_000, 0.1),
#         (100_000, 0.00001)
#     ]
    
#     for n, density in configs:
#         print(f"Generating matrix {n}x{n} with density {density}")
#         matrix = generate_sparse_coo_matrix(n, density)
#         filename = f"coo_input_data/matrix_{n}x{n}_d{int(density*100)}.mtx"
#         save_coo_matrix_market(matrix, filename)
#         print(f"Saved {filename} with {matrix.nnz} nonzeros")

# if __name__ == "__main__":
#     main()
