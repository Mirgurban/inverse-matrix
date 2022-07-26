#define GIGA_MODIFIER 1e9
#include <fstream>
#include <iostream>
#include <pthread.h>
#include <string.h>
#include <time.h>
#define args_4 5
#define args_5 6
#include "OBR.h"
#include "some_threads.h"
#include "vvod.h"
using namespace std;

struct thread_args {
    const float *a;
    const float *b;
    float *c;
    int n;
    int line;
};

int main(int argc, char *argv[])
{
    if (argc == 1) {
        cout << "Example of input" << endl;
        cout << "./a.out 20 4 4" << endl;
        return -1;
    }
    if (argc != args_4 && argc != args_5) {
        cout << "Wrong number of arguments" << endl;
        return -1;
    }
    int p = isnumber(argv[1]);
    int n = isnumber(argv[2]);
    int m = isnumber(argv[3]);
    int k = isnumber(argv[4]);
    int len_array = n * n;
    unsigned long long all_time = 0;

    float *array = NULL;

    if (n == -1 || m == -1 || k == -1) {
        cout << "cannot read number" << endl;
        return -1;
    }
    if (p < 1) {
        cout << "incorrect arguments" << endl;
        return -1;
    }
    if ((array = (float *)malloc(len_array * sizeof(float))) == nullptr) {
        cout << "failed to allocate memory" << endl;
        return -2;
    }
    unsigned long long *time_thread = new unsigned long long[p];
    for (int i = 0; i < p; ++i) {
        time_thread[i] = 0;
    }

    // ввод матрицы из файла или заполенение с помощью функции
    int er = vvod(array, n, m, k, argc, argv[args_4]);
    if (er != 0) {
        free(array);
        delete[] time_thread;
        return er;
    }

    vivod(array, n, n, m);
    cout << endl;

    // единичная матрица
    float *E = new float[n * n];
    fill_E(E, n);

    // изначальная матрица(нужна для невязки)
    float *old_array = new float[n * n];
    fill_old_array(old_array, array, n);

    all_time = currentTimeNano();
    // проверка на невырожденность матрицы
    // if (check_for_degeneration(n, array) == -4) {
    //     delete[] time_thread;
    //     delete[] old_array;
    //     delete[] E;
    //     return -4;
    // }

    // нахождение невязки
    float *mltp_mtrx = new float[n * n];
    fill_zeros(mltp_mtrx, n);
    float sum = thread_all(n, p, array, E, time_thread, mltp_mtrx, old_array);

    if (srav(sum, 0) == 1) {
        free(array);
        delete[] E;
        delete[] mltp_mtrx;
        delete[] old_array;
        delete[] time_thread;
        return -4;
    }
    cout << "Inverse: " << endl;
    vivod(E, n, n, m);
    cout << endl;

    all_time = currentTimeNano() - all_time;
    cout << "Residual: ";
    printf("%10.3e", sum);

    for (int i = 0; i < p; ++i) {
        cout << endl << "Time of thread " << i + 1 << ":";
        printf("%.2f", float(time_thread[i]) / GIGA_MODIFIER);
    }

    cout << endl << "Time:";
    printf("%.2f", float(all_time) / GIGA_MODIFIER);
    cout << endl;

    free(array);

    delete[] E;
    delete[] mltp_mtrx;
    delete[] old_array;
    delete[] time_thread;
    return 0;
}
