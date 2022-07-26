#include "math.h"
#include <iostream>
using namespace std;
#define eps 1.e-100
float *find_x_by_y(const float *y, int n)
{
    float pred_norma_y = 0;
    for (int i = 1; i < n; i++) {
        pred_norma_y += y[i] * y[i];
    }
    float module_y = sqrt(y[0] * y[0] + pred_norma_y);
    // module x'
    float x1 = y[0] - module_y;
    float module_x = sqrt(x1 * x1 + pred_norma_y);
    float *x = new float[n];
    x[0] = x1 / module_x;
    for (int i = 1; i < n; i++) {
        x[i] = y[i] / module_x;
    }
    return x;
}

float *U_times_y(const float *x, const float *y, int n)
{
    float scalar = 0;
    float *result = new float[n];
    for (int i = 0; i < n; i++) {
        scalar += y[i] * x[i];
    }
    for (int i = 0; i < n; i++) {
        result[i] = y[i] - 2 * scalar * x[i];
    }
    return result;
}

float *U_times_A(float *x, float *A, int k, int n)
{
    float *result = 0;

    for (int j = k; j < n; j++) {
        float y[n - k];
        for (int i = k; i < n; i++) {
            y[i - k] = A[i * n + j];
        }
        result = U_times_y(x, y, n - k);
        for (int i = k; i < n; i++) {
            A[i * n + j] = result[i - k];
        }
        delete[] result;
    }
    return A;
}

float *U_times_E(float *x, float *E, int k, int n)
{
    float *result = 0;

    for (int j = 0; j < n; j++) {
        float y[n - k];
        for (int i = k; i < n; i++) {
            y[i - k] = E[i * n + j];
        }
        result = U_times_y(x, y, n - k);
        for (int i = k; i < n; i++) {
            E[i * n + j] = result[i - k];
        }
        delete[] result;
    }
    return E;
}

float *multipl_mtrx(const float *a, const float *b, int n)
{
    float *bt = new float[n * n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            bt[i * n + j] = b[j * n + i];
        }
    }
    float *c = new float[n * n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            c[i * n + j] = 0.0;
            for (int k = 0; k < n; k++) {
                c[i * n + j] += a[i * n + k] * bt[j * n + k];
            }
        }
    }
    delete[] bt;
    return c;
}

float *gaus(float *E, const float *A, int n)
{
    for (int i = n - 1; i >= 0; i--) {
        float a = A[i * n + i];
        for (int j = 0; j < n; j++) {
            E[i * n + j] = E[i * n + j] / a;
            for (int l = 0; l < i; l++) {
                E[l * n + j] = E[l * n + j] - E[i * n + j] * A[l * n + i];
            }
        }
    }
    return E;
}

int srav(float a, float b)
{
    if (fabs(a - b) < eps) {
        return 1;
    }
    return 0;
}

void GetMatr(const float *A, float *p, int i, int n)
{
    int di = 0;
    int dj = 1;
    for (int ki = 0; ki < n - 1; ki++) {
        if (ki == i) {
            di = 1;
        }
        for (int kj = 0; kj < n - 1; kj++) {
            p[ki * (n - 1) + kj] = A[(ki + di) * n + kj + dj];
        }
    }
}

float Determinant(float *A, int n)
{
    float d = 0;
    float k = 1;
    float *p = new float[(n - 1) * (n - 1)];
    if (n == 1) {
        d = A[0];
        delete[] p;
        return d;
    }
    if (n == 2) {
        d = A[0] * A[3] - A[1] * A[2];
        delete[] p;
        return d;
    }
    if (n > 2) {
        for (int i = 0; i < n; i++) {
            GetMatr(A, p, i, n);
            d = d + k * A[i * n] * Determinant(p, n - 1);
            k = -k;
        }
    }
    delete[] p;
    return d;
}

float Residual(const float *A, int n)
{
    float sum = 0;
    for (int i = 0; i < n * n; ++i) {
        sum += A[i] * A[i];
    }
    sum = sqrt(sum);
    return sum;
}

void inverse_matrix(float *array, float *E, int n)
{
    for (int j = 0; j < n; j++) {
        int dimen = n - j;
        float y[dimen];
        int flag2 = 1;
        for (int i = 0; i < dimen; i++) {
            y[i] = array[n * (i + j) + j];
            if (i != 0 && srav(y[i], 0.0) == 0) {
                flag2 = 0;
            }
        }

        if (flag2 == 0) {
            float *x = find_x_by_y(y, dimen);
            array = U_times_A(x, array, j, n);
            E = U_times_E(x, E, j, n);
            delete[] x;
        }
    }
    gaus(E, array, n);
}

void fill_old_array(float *old_array, const float *array, int n)
{
    for (int i = 0; i < n * n; i++) {
        old_array[i] = array[i];
    }
}
