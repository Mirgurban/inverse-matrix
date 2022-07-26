#ifndef SOME_THREADS_H
#define SOME_THREADS_H

void *multipl_mtrx_line(void * args);
void multipl_mtrx(int p, const float *old_array, const float *E, float *mltp_mtrx, int n, unsigned long long *time_thread);
float *U_times_A(int p, const float *x, float *A, int k, int n);
float *U_times_E(int p, const float *x, float *A, int k, int n);
void *U_times_y_for_th_A(void * arg);
void *U_times_y_for_th_E(void * arg);
void inverse_matrix_th(int p, float *array, float *E, int n, unsigned long long *time_thread);
void *thread_all_help(void *arg);
float thread_all(int n, int p, float *A, float *E, unsigned long long *time_thread, float *mltp_mtrx, const float *old_array);
struct thread_args;

#endif  
