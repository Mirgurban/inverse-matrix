#define GIGA_MODIFIER 1e9
#include "math.h"
#include "vvod.h"
#include <iostream>
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <string.h>
#include <unistd.h>
using namespace std;
#define eps 1.e-15

// struct {
//     double val;
//     // int rank;
// } all_sum, sum;

unsigned long long currentTimeNano()
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (long long)(t.tv_sec * GIGA_MODIFIER + t.tv_nsec);
}

int srav(double a, double b)
{
    if (fabs(a - b) < eps) {
        return 1;
    }
    return 0;
}

double sum_2(const double *matrix_buffer, int n, int size)
{
    double sum = 0;
    for (int i = 0; i < n * (n / size + 1); ++i) {
        sum += matrix_buffer[i] * matrix_buffer[i];
    }
    return sum;
}

double *find_x_by_y(const double *y, int n, double *x)
{
    double pred_norma_y = 0;
    for (int i = 1; i < n; i++) {
        pred_norma_y += y[i] * y[i];
    }
    double module_y = sqrt(y[0] * y[0] + pred_norma_y);
    // module x'
    double x1 = y[0] - module_y;
    double module_x = sqrt(x1 * x1 + pred_norma_y);
    x[0] = x1 / module_x;
    for (int i = 1; i < n; i++) {
        x[i] = y[i] / module_x;
    }
    return x;
}

double *U_times_y(const double *x, const double *y, int n, double *result)
{
    double scalar = 0;
    for (int i = 0; i < n; i++) {
        scalar += y[i] * x[i];
    }
    for (int i = 0; i < n; i++) {
        result[i] = y[i] - 2 * scalar * x[i];
    }
    return result;
}

void *multipl_mtrx_line(int size, const double *matrix_buffer, const double *E,
                        double *matrix_c, int n, int rank)
{
    double *col_buffer = (double *)calloc(n / size + 1, sizeof(double));
    for (int col = 0; col < n; col++) {
        double *temp = (double *)calloc(n, sizeof(double));
        for (int row = 0; row < n; row += size) {
            MPI_Scatter(E + row + col / size * n, 1, MPI_DOUBLE,
                        col_buffer + row / size, 1, MPI_DOUBLE, col % size,
                        MPI_COMM_WORLD);
        }
        for (int row = 0; row < n; ++row) {
            for (int j = 0; j < n / size + 1; ++j) {
                temp[row] += matrix_buffer[row + j * n] * col_buffer[j];
            }
        }
        MPI_Reduce(temp, matrix_c + col / size * n, n, MPI_DOUBLE, MPI_SUM,
                   col % size, MPI_COMM_WORLD);
        delete[] temp;
    }

    matrix_minus_E(matrix_c, n, size, rank);

    delete[] col_buffer;
    return 0;
}

double Residual(const double *old_buffer, const double *E_buffer, int n,
                int size, int rank)
{
    double *matrix_c = (double *)calloc(n * (n / size + 1), sizeof(double));
    multipl_mtrx_line(size, old_buffer, E_buffer, matrix_c, n, rank);
    double sum = sum_2(matrix_c, n, size);
    double all_sum = 0;
    MPI_Reduce(&sum, &all_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        return sqrt(all_sum);
    }
    free(matrix_c);
    return 0;
}

double *U_times_A_column(double *x, double *matrix_buffer, int k, int n,
                         int column, int size)
{
    double *result = new double[n - k];
    double y[n - k];
    for (int i = k; i < n; i++) {
        y[i - k] = matrix_buffer[i + column / size * n];
    }
    U_times_y(x, y, n - k, result);
    for (int i = k; i < n; i++) {
        matrix_buffer[i + column / size * n] = result[i - k];
    }
    delete[] result;
    return matrix_buffer;
}

void GaussReverse(const double *matrix_column, double *E_buffer, int col, int n,
                  int size, int k)
{
    E_buffer[(n - k - 1) + col / size * n] =
        E_buffer[(n - k - 1) + col / size * n] / matrix_column[n - k - 1];
    double koeff = E_buffer[(n - k - 1) + col / size * n];
    for (int i = 0; i < n - k - 1; i++) {
        E_buffer[i + col / size * n] =
            E_buffer[i + col / size * n] - koeff * matrix_column[i];
    }
}

void *inverse_matrix(double *matrix_buffer, double *E_buffer, int n, int rank,
                     int size)
{
    int j0 = 0;
    for (int i = 0; i < n - 1; ++i) {
        j0 = i;
        double y[n - i];
        double *x = new double[n - i];
        int i_size_n = i / size * n;
        if (i % size == rank) {
            for (int h = i; h < n; h++) {
                y[h - i] = matrix_buffer[h + i_size_n];
            }
            find_x_by_y(y, n - i, x);
        }
        MPI_Bcast(x, n - i, MPI_DOUBLE, i % size, MPI_COMM_WORLD);
        for (int j = j0; j < n; ++j) {
            if (j % size == rank) {
                U_times_A_column(x, matrix_buffer, i, n, j, size);
            }
        }
        for (int j = 0; j < n; ++j) {
            if (j % size == rank) {
                U_times_A_column(x, E_buffer, i, n, j, size);
            }
        }
        delete[] x;
    }

    for (int k = 0; k < n; k++) {
        int size_n = size - n % size;
        double *matrix_column = new double[n - k]; // maybe not right
        // double matrix_column[n - k];
        if (rank == (n - k - 1) % size) {
            for (int i = 0; i < n - k; ++i) {
                matrix_column[i] = matrix_buffer[i + (n - k - 1) / size * n];
            }
        }
        MPI_Bcast(matrix_column, n - k, MPI_DOUBLE, (n - k - 1) % size,
                  MPI_COMM_WORLD);
        for (int col = 0; col < n + size_n; ++col) {
            if (col % size == rank) {
                GaussReverse(matrix_column, E_buffer, col, n, size, k);
            }
        }
        delete[] matrix_column;
    }
    return 0;
}