

// Compressed Sparse Row (CSR)
typedef struct {
    int n;             
    int nnz;  // Number of non-zero elements
    double* values;    // list of the non-zero vals
    int* col_index;    // Column index for each val
    int* row_ptr;      // 
} CSRMatrix;

// COO
typedef struct {
    int nnz;        // Number of non-zero elements
    int rows;       // Number of rows in the matrix
    int cols;       // Number of columns in the matrix
    int *row;       // Row indices (size nnz)
    int *col;       // Column indices (size nnz)
    double *values;    // Non-zero values (size nnz)
} COOMatrix;

COOMatrix create_coo(int nnz, int rows, int cols);

void free_coo(COOMatrix *mat);

COOMatrix read_and_create_coo(const char *filename);

void coo_matvec_mult(const COOMatrix *mat, const double *x, double *y);

void coo_matvec_mult_par(const COOMatrix *mat, const double *x, double *x_new, int start_idx, int local_nnz);

    void normalize_vector(double *x, int n);

void print_coo(const COOMatrix *mat);

void power_method_seq(const COOMatrix *mat, double *x, int max_iter, double tolerance);

void power_method_par(const COOMatrix *mat, double *x, int max_iter, double tolerance, int start_idx, int local_nnz);

void print_vector(double *vector, int n);

double dot_product(const double *a, const double *b, int n);

double compute_eigenvalue(const COOMatrix *mat, const double *x);

double validate_eigenpair(const COOMatrix *mat, const double *x, double lambda);

double vector_diff_norm(const double *a, const double *b, int n);