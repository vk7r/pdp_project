#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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

void power_method_seq(const COOMatrix *mat, double *x, int max_iter)
{
    int n = mat->rows;
    double *x_new = malloc(n * sizeof(double));



    for (int iter = 0; iter < max_iter; iter++)
    {

        // Ax
        coo_matvec_mult(mat, x, x_new);
        
        // x / ||x||
        normalize_vector(x_new, n);

        // copy x_new -> x
        for (int i = 0; i < n; i++)
        {
            // printf("x[%d] = %f\n", i, x[i]);
            x[i] = x_new[i];
        }

        
        // memcpy(x, x_new, n * sizeof(double));

        // print_vector(x_new, n);

    }

    // printf("x_new\n");
    // print_vector(x_new, n);

    // printf("x\n");
    // print_vector(x, n);

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