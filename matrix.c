#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "matrix.h"


COOMatrix create_coo(int nnz, int rows, int cols) {
    COOMatrix mat;
    mat.nnz = nnz;
    mat.rows = rows;
    mat.cols = cols;
    mat.row = (int *)malloc(nnz * sizeof(int));
    mat.col = (int *)malloc(nnz * sizeof(int));
    mat.values = (double *)malloc(nnz * sizeof(double));
    return mat;
}

void free_coo(COOMatrix *mat)
{
    free(mat->row);
    free(mat->col);
    free(mat->values);
}

COOMatrix read_and_create_coo(const char *filename) {
    FILE *f = fopen(filename, "r");

    // Skips first header line in file (%%MatrixMarket)
    char buffer[100];
    if (!fgets(buffer, sizeof(buffer), f))
    {
        fprintf(stderr, "Failed to read header line\n");
        fclose(f);
        exit(EXIT_FAILURE);
    }

    int rows, cols, nnz;

    if(fscanf(f, "%d %d %d", &rows, &cols, &nnz) != 3) {
        fprintf(stderr, "Error reading matrix dimensions\n");
        fclose(f);
        exit(EXIT_FAILURE);
    }

    COOMatrix mat = create_coo(nnz, rows, cols);

    // Populate the COO matrix
    for (int i = 0; i < nnz; ++i) {
        if (fscanf(f, "%d %d %lf", &mat.row[i], &mat.col[i], &mat.values[i]) != 3) {

            fprintf(stderr, "Error reading matrix data\n");
            free_coo(&mat);
            fclose(f);
            exit(EXIT_FAILURE);
        }

        // Convert to 0-based indexing
        mat.row[i]--;
        mat.col[i]--;
    }

    fclose(f);

    return mat;
}


void print_coo(const COOMatrix *mat) {
    printf("COO Matrix %dx%d:\n", mat->rows, mat->cols);
    for (int i = 0; i < mat->nnz; ++i) {
        printf("Row: %d, Col: %d, Value: %.2f\n", mat->row[i], mat->col[i], mat->values[i]);
    }
}

void coo_matvec_mult(const COOMatrix *mat, const double *x, double *x_new)
{
    // zero out the vector x_new
    for (int i = 0; i < mat->rows; i++)
    {
        x_new[i] = 0.0;
    }

    for (int i = 0; i < mat->nnz; i++)
    {
        x_new[mat->row[i]] += mat->values[i] * x[mat->col[i]];
    }

}

void coo_matvec_mult_par(const COOMatrix *mat, const double *x, double *x_new_local, int start_idx, int local_nnz)
{
    // for(int i = 0; i < local_nnz; i++)
    for (int i = start_idx; i < start_idx + local_nnz; i++)
    {
        int r = mat->row[i];
        int c = mat->col[i];
        double val = mat->values[i];

        x_new_local[r] += val * x[c];
    }
}

void normalize_vector(double *x, int n)
{
    double norm = 0.0;

    // caclulate the norm == sqrt(x[0]^2 + x[1]^2 + ... + x[n]^2)
    for (int i = 0; i < n; i++)
    {
        norm += x[i] * x[i];
    }
    
    norm = sqrt(norm);

    if (norm < 1e-10)
    {
        fprintf(stderr, "WARNING: Norm is too small, stopping power method early.\n");
        return;
    }

    // normalize the vector
    for (int i = 0; i < n; i++)
    {
        x[i] /= norm;
    }
}

// Helper function to compute difference norm between two vectors
double vector_diff_norm(const double *a, const double *b, int n)
{
    double sum = 0.0;
    for (int i = 0; i < n; i++)
    {
        double diff = a[i] - b[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}

void power_method_seq(const COOMatrix *mat, double *x, int max_iter, double tolerance)
{
    int n = mat->rows;
    double *x_new = malloc(n * sizeof(double));

    int iter = 0;
    double diff = tolerance + 1.0; // start bigger than tolerance to enter loop
    // printf("tolerance: %f\n", tolerance);
    int result = diff > tolerance;

    while (iter < max_iter && diff > tolerance)
    {

        // Ax
        coo_matvec_mult(mat, x, x_new);
        
        // x / ||x||
        normalize_vector(x_new, n);

        // Calculate difference norm between x_new and x_old
        diff = vector_diff_norm(x_new, x, n);
        // printf("Iteration %d, new diff: %e\n", iter, diff);

        // copy x_new -> x
        for (int i = 0; i < n; i++)
        {
            // printf("x[%d] = %f\n", i, x[i]);
            x[i] = x_new[i];
        }
        iter++;

    }

    free(x_new);
}

void power_method_par(const COOMatrix *mat, double *x, int max_iter, double tolerance, int start_idx, int local_nnz)
{
    int n = mat->rows;
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double *x_new_local = calloc(n, sizeof(double)); // partial Ax for this rank
    double *x_new = malloc(n * sizeof(double));      // combined Ax after allreduce

    int iter = 0;
    double diff = tolerance + 1.0;

    while (iter < max_iter && diff > tolerance)
    {


        // zero out local vector
        for (int i = 0; i < n; i++)
        {
            x_new_local[i] = 0.0;
        }

        // printf("Rank %d, IM DOING from %d to %d\n", rank, start_idx, start_idx + local_nnz);

        coo_matvec_mult_par(mat, x, x_new_local, start_idx, local_nnz);

        // Combine the local vector results of all ranks 
        MPI_Allreduce(x_new_local, x_new, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        normalize_vector(x_new, n);

        diff = vector_diff_norm(x_new, x, n);

        for (int i = 0; i < n; i++)
        {
            x[i] = x_new[i];
        }

        iter++;
    }

    free(x_new_local);
    free(x_new);
}


void print_vector(double *vector, int n)
{
    for (int i = 0; i < n; i++)
    {
        printf("%f ", vector[i]);
    }
    printf("\n");
}


double dot_product(const double *a, const double *b, int n)
{
    double sum = 0.0;
    for (int i = 0; i < n; i++)
    {
        sum += a[i] * b[i];
    }
    return sum;
}

// Compute eigenvalue lambda = (x^T A x) / (x^T x)
double compute_eigenvalue(const COOMatrix *mat, const double *x)
{
    int n = mat->rows;
    double *Ax = malloc(n * sizeof(double));

    coo_matvec_mult(mat, x, Ax);

    double numerator = dot_product(x, Ax, n);
    double denominator = dot_product(x, x, n);

    free(Ax);

    return numerator / denominator;
}

// Checks if the eigenpair (x, λ) is valid
double validate_eigenpair(const COOMatrix *mat, const double *x, double lambda)
{
    int n = mat->rows;
    double *Ax = malloc(n * sizeof(double));
    double residual = 0.0;

    coo_matvec_mult(mat, x, Ax); // Compute Ax

    for (int i = 0; i < n; i++)
    {
        double diff = Ax[i] - lambda * x[i]; // r = Ax - λx
        residual += diff * diff;
    }

    free(Ax);

    // printf("Checking residual norm ||Ax - λx||: %e\n", residual);
    if (residual < 1e-6)
    {
        printf("The eigenpair is CORRECT.\n");
    }
    else
    {
        printf("The eigenpair is NOT accurate!\n");

        exit(EXIT_FAILURE);
    }

    return sqrt(residual); // ||r||
}

