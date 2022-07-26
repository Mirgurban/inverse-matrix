#ifndef VVOD_H
#define VVOD_H

int vvod(float *array, int n, int m, int k, int argc, char *name);
unsigned long long currentTimeNano();
void vivod(float array[], int line, int column, int m);
int isnumber(char *name);
int error_check(int n, int m, int k, int num_of_arg);
float fill_matrix(int k, int n, int i, int j);
void fill_E(float * E, int n);

#endif  
