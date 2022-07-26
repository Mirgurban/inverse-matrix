#ifndef OBR_H
#define OBR_H

float* find_x_by_y(const float *y, int n, float *x);
float* U_times_y(const float *x, const float *y, int n, float *result);
int srav(float a,float b);
void fill_old_array(float *old_array, const float *array, int n);
int check_for_degeneration(int n, float *array);
struct thread_args;

#endif
