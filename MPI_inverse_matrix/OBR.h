#ifndef OBR_H
#define OBR_H

void *multipl_mtrx_line(int size, const double *matrix_buffer, const double *E, double *matrix_c,
                        int n, int rank);
int srav(double a, double b);
double *U_times_y(const double *x, const double *y, int n, double *result);
double *U_times_A_column(double *x, double *A, int k, int n, int column, int size);
double *find_x_by_y(const double *y, int n, double *x);
void *inverse_matrix(double *matrix_buffer, double *E_buffer, int n, int rank, int size);
// double Residual(const double *matrix_buffer, int n, int size, int rank);
double Residual(const double *old_buffer, const double *E_buffer, int n, int size, int rank);
unsigned long long currentTimeNano();

#endif  