#ifndef VVOD_H
#define VVOD_H

int isnumber(char *name);
double fill_matrix(int k, int n, int i, int j);
int read_row(int n, int m, int k, int argc, int row, double *row_buffer, char *name, int size);
void vivod(double *matrix_buffer, int size, int rank, int m, int n);
void fill_E(double *matrix_buffer, int n, int size, int rank);
void fill_old_array(double *old_array, const double *array, int n, int size);
void matrix_minus_E(double *matrix_buffer, int n, int size, int rank);

#endif  