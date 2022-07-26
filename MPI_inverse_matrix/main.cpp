#define GIGA_MODIFIER 1e9
#include "OBR.h"
#include "vvod.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <string.h>
#define args_4 4
#define args_5 5
using namespace std;

int main(int argc, char *argv[])
{
    if (argc == 1) {
        cout << "Example of input" << endl;
        cout << "mpirun -np 4 ./main 203 4 1" << endl;
        return -1;
    }
    if (argc != args_4 && argc != args_5) {
        cout << "Wrong number of arguments" << endl;
        return -1;
    }

    int n = isnumber(argv[1]);
    int m = isnumber(argv[2]);
    int k = isnumber(argv[3]);
    unsigned long long time = 0;
    unsigned long long time_all = 0;

    if (n == -1 || m == -1 || k == -1) {
        cout << "cannot read number" << endl;
        return -1;
    }

    if (n == 0 || m == 0) {
        cout << "incorrect arguments" << endl;
        return -1;
    }
    int rank;
    int size;
    time_all = currentTimeNano();
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double *matrix_buffer = NULL;
    double *row_buffer = NULL;

    if ((matrix_buffer = (double *)malloc(n * (n / size + 1) *
                                          sizeof(double))) == nullptr) {
        cout << "failed to allocate memory" << endl;
        return -2;
    }
    if ((row_buffer = (double *)malloc((n + size) * sizeof(double))) ==
        nullptr) {
        cout << "failed to allocate memory" << endl;
        free(matrix_buffer);
        return -2;
    }

    for (int row = 0; row < n; row++) {
        if (rank == 0) {
            int er = read_row(n, m, k, argc, row, row_buffer, argv[4], size);
            if (er != 0) {
                delete[] row_buffer;
                delete[] matrix_buffer;
                return er;
            }
        }
        for (int col = 0; col < n; col += size) {
            MPI_Scatter(row_buffer + col, 1, MPI_DOUBLE,
                        matrix_buffer + row + col / size * n, 1, MPI_DOUBLE, 0,
                        MPI_COMM_WORLD);
        }
    }

    vivod(matrix_buffer, size, rank, m, n);
    if (rank == 0) {
        cout << endl;
    }

    double *E_buffer = (double *)calloc(n * (n / size + 1), sizeof(double));
    double *old_buffer = (double *)calloc(n * (n / size + 1), sizeof(double));
    // double *matrix_c = (double *)calloc(n * (n / size + 1), sizeof(double));

    fill_E(E_buffer, n, size, rank);
    fill_old_array(old_buffer, matrix_buffer, n, size);

    if (rank == 0) {
        cout << "Inverse: " << endl;
    }
    // time = currentTimeNano();
    inverse_matrix(matrix_buffer, E_buffer, n, rank, size);
    // time = currentTimeNano() - time;
    for (int i = 0; i < (n / size + 1) * n; ++i) {
        if ((i % n == size * (i / n) + rank) &&
            (srav(matrix_buffer[i], 0) == 1)) {
            delete[] old_buffer;
            delete[] E_buffer;
            delete[] matrix_buffer;
            // delete[] matrix_c;
            delete[] row_buffer;
            cout << "degenerate matrix";
            return -4;
        }
    }

    // multipl_mtrx_line(size, old_buffer, E_buffer, matrix_c, n, rank);
    // matrix_minus_E(matrix_c, n, size, rank);
    vivod(E_buffer, size, rank, m, n);

    if (rank == 0) {
        cout << endl;
    }
    time = currentTimeNano();
    double Res = Residual(old_buffer, E_buffer, n, size, rank);
    time = currentTimeNano() - time;
    if (rank == 0) {
        cout << "Residual: ";
        printf("%10.3e\n", Res);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    printf("%s%d%s%.2f\n", "Time of process ", rank, ": ",
           float(time) / GIGA_MODIFIER);

    MPI_Finalize();
    time_all = currentTimeNano() - time_all;
    if (rank == 0) {
        printf("%s%.2f\n", "Time:", float(time_all) / GIGA_MODIFIER);
    }
    delete[] old_buffer;
    delete[] E_buffer;
    delete[] matrix_buffer;
    // delete[] matrix_c;
    delete[] row_buffer;

    return 0;
}
