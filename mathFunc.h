#include <cmath>
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>
#include <thread>
#include <chrono>
#include <dirent.h>
#include <random>
#include <iostream>
#include <algorithm>  
#include <string.h>
#include <iomanip>
#include "adios2.h"
#include "mgard/compress_x.hpp"


int idx(int i , int j , int N) {
    return i * N + j;
}

// X-direction central difference
std::vector<double> centralDifferenceX(const std::vector<double>& f , int N , double elementLength) {
    std::vector<double> diff(N * N , 0.0);

    for (int i = 1; i < N - 1; ++i) {
        for (int j = 1; j < N - 1; ++j) {
            int center = idx(i , j , N);
            int left = idx(i , j - 1 , N);
            int right = idx(i , j + 1 , N);

            diff[center] = (f[right] - f[left]) / (2.0 * elementLength);
        }
    }

    return diff;
}

// Y-direction central difference
std::vector<double> centralDifferenceY(const std::vector<double>& f , int N , double elementLength) {
    std::vector<double> diff(N * N , 0.0);

    for (int i = 1; i < N - 1; ++i) {
        for (int j = 1; j < N - 1; ++j) {
            int center = idx(i , j , N);
            int up = idx(i - 1 , j , N);
            int down = idx(i + 1 , j , N);

            diff[center] = (f[down] - f[up]) / (2.0 * elementLength);
        }
    }

    return diff;
}


// --------------------------------------------------------------------------

// Laplaccian method
std::vector<double> laplacian(
    const std::vector<double>& f ,
    int N ,
    double elementLength
) {
    std::vector<double> laplace(N * N , 0.0);

    for (int i = 1; i < N - 1; ++i) {
        for (int j = 1; j < N - 1; ++j) {
            int center = idx(i , j , N);
            int up = idx(i - 1 , j , N);
            int down = idx(i + 1 , j , N);
            int left = idx(i , j - 1 , N);
            int right = idx(i , j + 1 , N);

            laplace[center] = (
                f[up] + f[down] + f[left] + f[right] - 4.0 * f[center]
                ) / (elementLength * elementLength);
        }
    }

    return laplace;
}