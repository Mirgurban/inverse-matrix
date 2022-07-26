#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#define GIGA_MODIFIER 1e9
#include "math.h"
#include <fstream>
#include <iostream>
#include <string.h>
#define args_4 4
#define args_5 5
using namespace std;

double fill_matrix(int k, int n, int i, int j)
{
    switch (k) {
    case 1: {
        return (double)(n - max(i + 1, j + 1) + 1);
    }
    case 2: {
        return (double)(max(i + 1, j + 1));
    }
    case 3: {
        return (double)(abs(i - j));
    }
    case 4: {
        return (double)(1.0 / (i + j + 1));
    }
    }
    return -1;
}

int isnumber(char *name)
{
    for (long unsigned int i = 0; i < strlen(name); ++i) {
        int flag = 0;
        if (isdigit(name[i]) == 0) {
            flag = 1;
        }
        if (flag == 1) {
            return -1;
        }
    }
    return atoi(name);
}

int error_check(int n, int m, int k, int num_of_arg)
{
    if (n <= 0) {
        cout << "incorrect arguments" << endl;
        return -1;
    }
    if (m <= 0) {
        cout << "incorrect arguments" << endl;
        return -1;
    }
    if (m > n) {
        cout << "incorrect arguments" << endl;
        return -1;
    }
    if (k < 0 || k > 4) {
        cout << "incorrect arguments" << endl;
        return -1;
    }
    if (num_of_arg == args_5 && k != 0) {
        cout << "incorrect arguments" << endl;
        return -1;
    }
    if (k == 0 && num_of_arg != args_5) {
        cout << "incorrect arguments" << endl;
        return -1;
    }
    return 0;
}

int read_row(int n, int m, int k, int argc, int row, double *row_buffer,
             char *name, int size)
{
    double number;
    if (error_check(n, m, k, argc) != 0) {
        return -1;
    }
    // check for argc and errors
    if (k != 0) {
        for (int i = 0; i < n; ++i) {
            row_buffer[i] = fill_matrix(k, n, row, i);
        }
        // for (int i = n; i < n + size; ++i) {
        //     row_buffer[i] = 0;
        // }
    } else {
        ifstream file;
        file.open(name);
        if (!(file.is_open())) {
            cout << "cannot open file" << endl;
            return -3;
        }
        for (int i = 0; i < n * row; ++i) {
            file >> number;
        }
        for (int i = 0; i < n; i++) {
            file >> number;
            row_buffer[i] = number;
            if (!file.good()) {
                cout << "something wrong with numbers in file" << endl;
                return -3;
            }
        }
        for (int i = n; i < n + size; ++i) {
            row_buffer[i] = 0;
        }
        file.close();
    }
    return 0;
}

void vivod(double *matrix_buffer, int size, int rank, int m, int n)
{
    double *row_buffer = (double *)calloc(n + size, sizeof(double));
    for (int row = 0; row < m; ++row) {
        for (int col = 0; col < m; col += size) {
            MPI_Gather(matrix_buffer + row + col / size * n, 1, MPI_DOUBLE,
                       row_buffer + col, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        if (rank == 0) {
            for (int i = 0; i < m; ++i) {
                printf("%10.3e", row_buffer[i]);
            }
            cout << endl;
        }
    }
    delete[] row_buffer;
}

void fill_E(double *matrix_buffer, int n, int size, int rank)
{
    int num_of_el = (n / size + 1) * n;
    for (int i = 0; i < num_of_el; ++i) {
        if (i % n == size * (i / n) + rank) {
            matrix_buffer[i] = 1.0;
        } else {
            matrix_buffer[i] = 0.0;
        }
    }
}

void fill_old_array(double *old_array, const double *array, int n, int size)
{
    for (int i = 0; i < n * (n / size + 1); i++) {
        old_array[i] = array[i];
    }
}

void matrix_minus_E(double *matrix_buffer, int n, int size, int rank)
{
    int num_of_el = (n / size + 1) * n;
    for (int i = 0; i < num_of_el; ++i) {
        if (i % n == size * (i / n) + rank) {
            matrix_buffer[i] = matrix_buffer[i] - 1.0;
        }
    }
}