#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/stat.h>

typedef struct {
    int nnz;
    int N;
    int *row;
    int *col;
    double *val;
} COOMatrix;

COOMatrix generate_random_coo(int N, double density) {
    int nnz = (int)(density * N * N);
    COOMatrix mat;
    mat.N = N;
    mat.nnz = nnz;
    mat.row = malloc(nnz * sizeof(int));
    mat.col = malloc(nnz * sizeof(int));
    mat.val = malloc(nnz * sizeof(double));

    for (int i = 0; i < nnz; ++i) {
        mat.row[i] = rand() % N;
        mat.col[i] = rand() % N;
        mat.val[i] = ((double)rand() / RAND_MAX) * 200.0 - 100.0; // random in [-1, 1]
    }

    return mat;
}

void write_matrix_market(const COOMatrix *mat, const char *filename) {
    FILE *f = fopen(filename, "w");
    if (!f) {
        perror("File open failed");
        exit(EXIT_FAILURE);
    }

    fprintf(f, "%%MatrixMarket matrix coordinate real general\n");
    fprintf(f, "%d %d %d\n", mat->N, mat->N, mat->nnz);

    for (int i = 0; i < mat->nnz; ++i) {
        fprintf(f, "%d %d %.10f\n", mat->row[i] + 1, mat->col[i] + 1, mat->val[i]); // 1-based
    }

    fclose(f);
}

void free_coo(COOMatrix *mat) {
    free(mat->row);
    free(mat->col);
    free(mat->val);
}

int main() {
    srand(time(NULL));

    // Create output folder if not exists
    struct stat st = {0};
    if (stat("coo_input_data", &st) == -1) {
        mkdir("coo_input_data", 0700);
    }

    // Define test sizes and densities
    int sizes[] = {50, 100, 1000, 5000, 10000};
    // int sizes[] = {10, 50, 100};
    double densities[] = {0.0005, 0.001, 0.005};  // 0.05%, 0.1%, 0.5%

    for (int i = 0; i < sizeof(sizes)/sizeof(sizes[0]); ++i) {
        for (int j = 0; j < sizeof(densities)/sizeof(densities[0]); ++j) {
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
