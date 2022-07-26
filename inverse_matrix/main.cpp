#define GIGA_MODIFIER 1e9
#include <fstream>
#include <iostream>
#include <string.h>
#include <time.h>
#define args_4 4
#define args_5 5
#include "OBR.h"
#include "vvod.h"
using namespace std;

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

    int n = isnumber(argv[1]);
    int m = isnumber(argv[2]);
    int k = isnumber(argv[3]);
    int len_array = n * n;
    unsigned long long time = 0;
    float *array = NULL;

    if (n == -1 || m == -1 || k == -1) {
        cout << "cannot read number" << endl;
        return -1;
    }

    if ((array = (float *)malloc(len_array * sizeof(float))) == nullptr) {
        cout << "failed to allocate memory" << endl;
        return -2;
    }

    //ввод матрицы из файла или заполенение с помощью функции
    int er = vvod(array, n, m, k, argc, argv[4]);
    if (er != 0) {
        free(array);
        return er;
    }

    vivod(array, n, n, m);
    cout << endl;

    //единичная матрица
    float *E = new float[n * n];
    fill_E(E, n);

    //изначальная матрица(нужна для невязки)
    float *old_array = new float[n * n];
    fill_old_array(old_array, array, n);

    //в Е записывает обратную матрицу
    time = currentTimeNano();
    inverse_matrix(array, E, n);

    vivod(array, n,n,m);
    //проверка на невырожденность матрицы
    for (int i = 0; i < n; i++) {
        if (srav(array[i * n + i], 0) == 1) {
            free(array);
            delete[] E;
            delete[] old_array;
            cout << "degenerate matrix";
            return -4;
        }
    }

    cout << "Inverse: " << endl;
    vivod(E, n, n, m);
    cout << endl;
    //нахождение невязки
    float *mtrx = 0;
    mtrx = multipl_mtrx(old_array, E, n);

    for (int i = 0; i < n; i++) {
        mtrx[i * n + i] = mtrx[i * n + i] - 1;
    }

    cout << "Residual: ";
    printf("%10.3e", Residual(mtrx, n));
    time = currentTimeNano() - time;
    cout << endl << "Time:";
    printf("%.2f", float(time) / GIGA_MODIFIER);

    free(array);
    delete[] E;
    delete[] mtrx;
    delete[] old_array;
    return 0;
}
