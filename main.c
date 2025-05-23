#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "matrix.h"
#include "power_method.h"

#define MAX_ITER 1000000000
#define TOLERANCE 1e-8

// CSR or COO format ? --> CSR the better

// Parses the input 
// sets the values of the variables
// Calls the Power method
// Print results

// Distribute from root function??

void distribute_from_root(COOMatrix *mat, int rank, int size)
{
    COOMatrix full_mat;
    int *sendcounts = malloc(size * sizeof(int));
    int *chunk_starts = malloc(size * sizeof(int));

    if (rank == 0)
    {
        full_mat = *mat;
        int chunk_size = full_mat.nnz / size;
        int remainder = full_mat.nnz % size;
        int offset = 0;
        for (int i = 0; i < size; i++)
        {
            sendcounts[i] = (i < remainder) ? chunk_size + 1 : chunk_size;
            chunk_starts[i] = offset;
            offset += sendcounts[i];
        }
    }

    int nrows, ncols;
    if (rank == 0)
    {
        nrows = full_mat.rows;
        ncols = full_mat.cols;
    }
    MPI_Bcast(&nrows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ncols, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int local_nnz;
    MPI_Scatter(sendcounts, 1, MPI_INT, &local_nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);

    mat->nnz = local_nnz;
    mat->rows = nrows;
    mat->cols = ncols;

    mat->row = malloc(local_nnz * sizeof(int));
    mat->col = malloc(local_nnz * sizeof(int));
    mat->values = malloc(local_nnz * sizeof(double));

    MPI_Scatterv(full_mat.row, sendcounts, chunk_starts, MPI_INT,
                 mat->row, local_nnz, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(full_mat.col, sendcounts, chunk_starts, MPI_INT,
                 mat->col, local_nnz, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(full_mat.values, sendcounts, chunk_starts, MPI_DOUBLE,
                 mat->values, local_nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(sendcounts);
}

int main(int argc, char **argv)
{

    MPI_Init(&argc, &argv);

    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc < 2)
    {
        if (rank == 0)
        {
            fprintf(stderr, "Usage: %s <input_matrix>\n", argv[0]);
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    const char *input_matrix_file = argv[1];

    // ===============================================================
    //                           SEQUENTIAL SOLUTION
    // ===============================================================
    if(size == 1)
    {
        // printf("Running sequential...\n\n");
        COOMatrix mat = read_and_create_coo(input_matrix_file);
        
        int n = mat.rows;
        int max_iterations = MAX_ITER;
        double tolerance = TOLERANCE;

        // Create the vector x, with 1.0 in all positions
        double *x = malloc(n * sizeof(double));
        int vector_size = n;
        
        for (int i = 0; i < n; i++)
        {
            x[i] = 1.0; //((double)rand() / RAND_MAX) * 200.0 - 100.0;
        }
        
        double start = MPI_Wtime();
        power_method_seq(&mat, x, max_iterations, tolerance);
        double execution_time = MPI_Wtime() - start;
        
        // printf("Finished power method.\n\n");

        // print_vector(x, vector_size);

        double eigenvalue = compute_eigenvalue(&mat, x);
        printf("Final dominant eigenvalue: %f\n\n", eigenvalue);

        // Validate eigenpair
        double res = validate_eigenpair(&mat, x, eigenvalue);

        printf("execution time: %f seconds\n\n", execution_time);

        free_coo(&mat);
        MPI_Finalize();
        return 0;
    }

    // ===============================================================
    //                           PARALLEL SOLUTION
    // ===============================================================
    else
    {

        COOMatrix mat;
        COOMatrix full_mat;

        if (rank == 0)
        {
            full_mat = read_and_create_coo(input_matrix_file);
            mat = full_mat; // mat gets overwritten later anyway
        }

        distribute_from_root(&mat, rank, size);

        // printf("Rank %d: chunk = [%d, %d] local_nnz = %d\n", rank, start_idx, start_idx + mat.nnz, mat.nnz);

        int local_nnz = mat.nnz;
        int n = mat.rows;

        double *x = malloc(n * sizeof(double));
        int vector_size = n;
        int max_iterations = MAX_ITER;
        double tolerance = TOLERANCE;

        for (int i = 0; i < n; i++)
        {
            x[i] = 1.0; //((double)rand() / RAND_MAX) * 200.0 - 100.0;
        }

        double start = MPI_Wtime();

        power_method_par(&mat, x, max_iterations, tolerance);

        double execution_time = MPI_Wtime() - start;

        double max_time;
        MPI_Reduce(&execution_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        if(rank == 0)
        {
            // printf("Finished power method.\n\n");

            // print_vector(x, vector_size);

            double eigenvalue = compute_eigenvalue(&full_mat, x);
            printf("Final Dominant eigenvalue: %f\n", eigenvalue);

            // Validate eigenpair
            double res = validate_eigenpair(&full_mat, x, eigenvalue);

            printf("\nMax execution time: %f seconds\n\n", execution_time);

            free_coo(&full_mat);
        }
        
        // print_coo(&mat);
        free_coo(&mat);
        MPI_Finalize();
        return 0;
    }
}
