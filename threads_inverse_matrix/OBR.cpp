#include "math.h"
#include <iostream>
using namespace std;
#define eps 1.e-100

struct thread_args {
    const float *a;
    const float *b;
    float *c;
    int n;
    int line;
};

float *find_x_by_y(const float *y, int n, float *x)
{
    float pred_norma_y = 0;
    for (int i = 1; i < n; i++) {
        pred_norma_y += y[i] * y[i];
    }
    float module_y = sqrt(y[0] * y[0] + pred_norma_y);
    // module x'
    float x1 = y[0] - module_y;
    float module_x = sqrt(x1 * x1 + pred_norma_y);
    // float *x = new float[n];
    x[0] = x1 / module_x;
    for (int i = 1; i < n; i++) {
        x[i] = y[i] / module_x;
    }
    return x;
}

float *U_times_y(const float *x, const float *y, int n, float *result)
{
    float scalar = 0;
    // float *result = new float[n];
    for (int i = 0; i < n; i++) {
        scalar += y[i] * x[i];
    }
    for (int i = 0; i < n; i++) {
        result[i] = y[i] - 2 * scalar * x[i];
    }
    return result;
}

int srav(float a, float b)
{
    if (fabs(a - b) < eps) {
        return 1;
    }
    return 0;
}

void fill_old_array(float *old_array, const float *array, int n)
{
    for (int i = 0; i < n * n; i++) {
        old_array[i] = array[i];
    }
}

int check_for_degeneration(int n, float *array)
{
    for (int i = 0; i < n; i++) {
        if (srav(array[i * n + i], 0) == 1) {
            // free(array);
            cout << "degenerate matrix";
            return -4;
        }
    }
    return 0;
}
