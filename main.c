#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "matrix.h"
#include "power_method.h"

// CSR or COO format ? --> CSR the better

// Parses the input 
// sets the values of the variables
// Calls the Power method
// Print results
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

    COOMatrix mat = read_and_create_coo(input_matrix_file);

    int n = mat.rows;
    int iterations = 10;

    // Create the vector x, with 1.0 in all positions 
    double *x = malloc(n * sizeof(double));
    int vector_size = n;
    
    for (int i = 0; i < n; i++)
    {
        x[i] = 1.0; //((double)rand() / RAND_MAX) * 200.0 - 100.0;
    }


    //________________________ SEQUENTIAL SOLUTION ________________________//
    //_____________________________________________________________________//
    if(size == 1)
    {
        printf("Running sequential...\n");
    
        // print_coo(&mat);
        
        // printf("Initial vector:\n");
        // print_vector(x, vector_size);

        printf("Performing power method for %d iterations:\n", iterations);
        power_method_seq(&mat,x,iterations);

        printf("power method finished.\n");
        // Print result vector
        // printf("Final vector:\n");
        print_vector(x, vector_size);


        free_coo(&mat);
        MPI_Finalize();
        return 0;
    }


    //________________________ PARALLEL SOLUTION __________________________//
    //_____________________________________________________________________//
    else
    {
        printf("Running parallel...\n");
        // Parallel code here
        print_coo(&mat);
        free_coo(&mat);
        MPI_Finalize();
        return 0;
    }
}
