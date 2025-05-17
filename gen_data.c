#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/stat.h>

typedef struct
{
    int nnz;
    int N;
    int *row;
    int *col;
    double *val;
} COOMatrix;

// Helper function to check for duplicates
int exists(int *row, int *col, int count, int r, int c)
{
    for (int i = 0; i < count; ++i)
    {
        if (row[i] == r && col[i] == c)
            return 1;
    }
    return 0;
}

COOMatrix generate_random_coo(int N, double density)
{
    int max_possible_nnz = N * N;
    int nnz = (int)(density * max_possible_nnz);
    if (nnz > max_possible_nnz)
        nnz = max_possible_nnz;

    COOMatrix mat;
    mat.N = N;
    mat.nnz = nnz;
    mat.row = malloc(nnz * sizeof(int));
    mat.col = malloc(nnz * sizeof(int));
    mat.val = malloc(nnz * sizeof(double));

    int count = 0;

    // Ensure every row has at least one entry
    for (int i = 0; i < N && count < nnz; ++i)
    {
        int r = i;
        int c = rand() % N;
        if (!exists(mat.row, mat.col, count, r, c))
        {
            mat.row[count] = r;
            mat.col[count] = c;
            mat.val[count] = ((double)rand() / RAND_MAX) * 2.0 - 1.0; // [-1, 1]
            ++count;
        }
    }

    // Fill the rest randomly avoiding duplicates
    while (count < nnz)
    {
        int r = rand() % N;
        int c = rand() % N;
        if (!exists(mat.row, mat.col, count, r, c))
        {
            mat.row[count] = r;
            mat.col[count] = c;
            mat.val[count] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
            ++count;
        }
    }

    return mat;
}

void write_matrix_market(const COOMatrix *mat, const char *filename)
{
    FILE *f = fopen(filename, "w");
    if (!f)
    {
        perror("File open failed");
        exit(EXIT_FAILURE);
    }

    fprintf(f, "%%MatrixMarket matrix coordinate real general\n");
    fprintf(f, "%d %d %d\n", mat->N, mat->N, mat->nnz);

    for (int i = 0; i < mat->nnz; ++i)
    {
        fprintf(f, "%d %d %.10f\n", mat->row[i] + 1, mat->col[i] + 1, mat->val[i]); // 1-based
    }

    fclose(f);
}

void free_coo(COOMatrix *mat)
{
    free(mat->row);
    free(mat->col);
    free(mat->val);
}

int main()
{
    srand(time(NULL));

    // Create output folder if it doesn't exist
    struct stat st = {0};
    if (stat("coo_input_data", &st) == -1)
    {
        mkdir("coo_input_data", 0700);
    }

    // Define test sizes and densities
    int sizes[] = {50, 100, 1000, 5000, 10000};
    double densities[] = {0.0005, 0.001, 0.005}; // 0.05%, 0.1%, 0.5%

    for (int i = 0; i < sizeof(sizes) / sizeof(sizes[0]); ++i)
    {
        for (int j = 0; j < sizeof(densities) / sizeof(densities[0]); ++j)
        {
            int N = sizes[i];
            double d = densities[j];

            printf("Generating %dx%d with density %.4f...\n", N, N, d);

            COOMatrix mat = generate_random_coo(N, d);

            char filename[256];
            snprintf(filename, sizeof(filename),
                     "coo_input_data/input_%dx%d_d%.4f.mtx", N, N, d);

            write_matrix_market(&mat, filename);
            free_coo(&mat);
        }
    }

    printf("Done.\n");
    return 0;
}
