#define GIGA_MODIFIER 1e9
#include <fstream>
#include <iostream>
#include <string.h>
#include <unistd.h>
#define args_4 4
#define args_5 5
using namespace std;

float fill_matrix(int k, int n, int i, int j)
{
    switch (k) {
    case 1: {
        return (float)(n - max(i + 1, j + 1) + 1);
    }
    case 2: {
        return (float)(max(i + 1, j + 1));
    }
    case 3: {
        return (float)(abs(i - j));
    }
    case 4: {
        return (float)(1.0 / (i + j + 1));
    }
    }
    return -1;
}

int error_check(int n, int m, int k, int num_of_arg)
{
    if (n <= 0) {
        cout << "incorrect arguments" << endl;
        return -1;
    }
    if (m <= 0) {
        cout << "incorrect arguments" << endl;
        return -1;
    }
    if (m > n) {
        cout << "incorrect arguments" << endl;
        return -1;
    }
    if (k < 0 || k > 4) {
        cout << "incorrect arguments" << endl;
        return -1;
    }
    if (num_of_arg == args_5 && k != 0) {
        cout << "incorrect arguments" << endl;
        return -1;
    }
    if (k == 0 && num_of_arg != args_5) {
        cout << "incorrect arguments" << endl;
        return -1;
    }
    return 0;
}

int vvod(float *array, int n, int m, int k, int argc, char *name)
{
    float number;
    if (error_check(n, m, k, argc) != 0) {
        return -1;
    }
    if (k == 0 && argc == args_5) {
        ifstream file;
        file.open(name);
        if (!(file.is_open())) {
            cout << "cannot open file" << endl;
            return -3;
        }
        for (int i = 0; i < n * n; i++) {
            file >> number;
            array[i] = number;
            if (!file.good()) {
                cout << "something wrong with numbers in file" << endl;
                return -3;
            }
        }
        file.close();
    } else {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                array[n * i + j] = fill_matrix(k, n, i, j);
            }
        }
    }
    return 0;
}

unsigned long long currentTimeNano()
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (long long)(t.tv_sec * GIGA_MODIFIER + t.tv_nsec);
}

void vivod(float array[], int line, int column, int m)
{
    int i;
    int j;
    int min_line = min(m, line);
    int min_column = min(m, column);
    for (i = 0; i < min_line; i++) {
        for (j = 0; j < min_column; j++) {
            printf("%10.3e", array[line * i + j]);
        }
        cout << "\r\n";
    }
}

int isnumber(char *name)
{
    for (long unsigned int i = 0; i < strlen(name); ++i) {
        int flag = 0;
        if (isdigit(name[i]) == 0) {
            flag = 1;
        }
        if (flag == 1) {
            return -1;
        }
    }
    return atoi(name);
}

void fill_E(float *E, int n)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                E[i * n + i] = 1.0;
            } else {
                E[i * n + j] = 0.0;
            }
        }
    }
}
