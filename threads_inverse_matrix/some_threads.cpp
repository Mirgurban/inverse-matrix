#define GIGA_MODIFIER 1e9
#include "OBR.h"
#include "math.h"
#include "vvod.h"
#include <iostream>
using namespace std;

struct thread_args {
    float *c;
    float *A;
    float *E;
    int n;
    int p;
    int num_of_thread;
    unsigned long long time_thread;
    const float *old_array;
    int flag;
    pthread_barrier_t *barrier;
};

unsigned long long currentTimeNanoThread()
{
    struct timespec t;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
    return (long long)(t.tv_sec * GIGA_MODIFIER + t.tv_nsec);
}

void *multipl_mtrx_line(int p, const float *old_array, const float *E, float *c,
                        int n, int num_of_thread)
{
    int i = 0;
    while ((num_of_thread + i * p) * n + n - 1 < n * n) {
        for (int j = 0; j < n; j++) {
            float s = 0.0;
            for (int k = 0; k < n; k++) {
                s += old_array[(num_of_thread + i * p) * n + k] * E[k * n + j];
            }
            c[i * n + j] = s;
        }
        ++i;
    }
    return 0;
}

void GaussReverse(const float *A, float *E, int col, int n)
{
    for (int k = 1; k <= n; k++) {
        float koeff = E[(n - k) * n + col];
        for (int i = n - k + 1; i < n; i++) {
            koeff -= E[i * n + col] * A[(n - k) * n + i];
        }
        E[(n - k) * n + col] = koeff / A[(n - k) * n + n - k];
    }
}

float *U_times_A_column(float *x, float *A, int k, int n, int column)
{
    float *result = new float[n - k];
    int j = column;
    float y[n - k];
    for (int i = k; i < n; i++) {
        y[i - k] = A[i * n + j];
    }
    U_times_y(x, y, n - k, result);
    for (int i = k; i < n; i++) {
        A[i * n + j] = result[i - k];
    }
    delete[] result;
    return A;
}

void *inverse_matrix_th_help(int p, float *A, float *E, int n,
                             int num_of_thread, pthread_barrier_t *barrier)
{
    int j0 = 0;
    for (int i = 0; i < n - 1; ++i) {
        j0 = i;
        float y[n - i];
        for (int h = i; h < n; h++) {
            y[h - i] = A[h * n + i];
        }
        float *x = new float[n - i];
        find_x_by_y(y, n - i, x);
        pthread_barrier_wait(barrier);
        for (int j = j0; j < n; ++j) {
            if (j % p == num_of_thread) {
                U_times_A_column(x, A, i, n, j);
            }
        }

        for (int j = 0; j < n; ++j) {
            if (j % p == num_of_thread) {
                U_times_A_column(x, E, i, n, j);
            }
        }
        pthread_barrier_wait(barrier);
        delete[] x;
    }
    for (int k = 0; k < n; k++) {
        if (k % p == num_of_thread) {
            GaussReverse(A, E, k, n);
        }
    }
    return 0;
}

void *thread_all_help(void *arg)
{
    thread_args *args = (thread_args *)arg;
    int flag = args->flag;
    float *c = args->c;
    int n = args->n;
    int p = args->p;
    unsigned long long time_thread = args->time_thread;
    const float *old_array = args->old_array;
    float *A = args->A;
    float *E = args->E;
    int num_of_thread = args->num_of_thread;
    pthread_barrier_t *barrier = args->barrier;

    time_thread = currentTimeNanoThread();
    inverse_matrix_th_help(p, A, E, n, num_of_thread, barrier);
    pthread_barrier_wait(barrier);
    // проверка на невырожденность матрицы
    // vivod(A, n, n ,n);
    if (check_for_degeneration(n, A) == -4) {
        // cout << "deg " << num_of_thread << endl;
        flag = 1;
    }
    if (flag == 0) {
        multipl_mtrx_line(p, old_array, E, c, n, num_of_thread);
        args->time_thread = currentTimeNanoThread() - time_thread;
    }
    args->flag = flag;
    return 0;
}

float thread_all(int n, int p, float *A, float *E,
                 unsigned long long *time_thread, float *mltp_mtrx,
                 const float *old_array)
{
    pthread_t *thread = new pthread_t[p];
    thread_args *args = new thread_args[p];

    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, p);

    for (int i = 0; i < p; ++i) {
        args[i].n = n;
        args[i].E = E;
        args[i].A = A;
        args[i].old_array = old_array;
        args[i].c = new float[n * (n / p + 1)];
        args[i].barrier = &barrier;
        args[i].p = p;
        args[i].num_of_thread = i;
        args[i].time_thread = 0;
        args[i].flag = 0;
    }
    for (int i = 0; i < p; ++i) {
        pthread_create(&thread[i], NULL, thread_all_help, &args[i]);
    }
    for (int i = 0; i < p; ++i) {
        pthread_join(thread[i], NULL);
        time_thread[i] = args[i].time_thread;
    }
    pthread_barrier_destroy(&barrier);
    float sum = 0;
    for (int i = 0; i < p; ++i) {
        // cout << args[i].flag << endl;
        if (args[i].flag == 1) {
            delete[] thread;
            for (int j = 0; j < p; ++j) {
                delete[] args[j].c;
            }
            delete[] args;
            return sum;
        }
    }

    for (int i = 0; i < p; ++i) {
        for (int k = 0; k < n / p + 1; ++k) {
            for (int j = 0; j < n; ++j) {
                if ((i + p * k) * n + j < n * n) {
                    mltp_mtrx[(i + p * k) * n + j] = args[i].c[k * n + j];
                }
            }
        }
    }

    for (int i = 0; i < n; i++) {
        mltp_mtrx[i * n + i] = mltp_mtrx[i * n + i] - 1;
    }
    for (int i = 0; i < n * n; ++i) {
        sum += mltp_mtrx[i] * mltp_mtrx[i];
    }
    sum = sqrt(sum);

    delete[] thread;
    for (int i = 0; i < p; ++i) {
        delete[] args[i].c;
    }
    delete[] args;
    return sum;
}
