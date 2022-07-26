#ifndef OBR_H
#define OBR_H

float* find_x_by_y(const float *y, int n);
float* U_times_y(const float *x, const float *y, int n);
float* U_times_A(float x[], float *A, int k, int n);
float * multipl_mtrx(const float *a, const float *b, int n);
float* gaus(float *E, const float *A, int n);
int srav(float a,float b);
float* U_times_E(float *x, float *E, int k, int n);
float Determinant(float *A, int n);
void GetMatr(float *A, float *p, int i, int j, int n);
void inverse_matrix(float *array, float *E, int n);
void fill_old_array(float *old_array, const float *array, int n);
float Residual(const float *A, int n);

#endif  
