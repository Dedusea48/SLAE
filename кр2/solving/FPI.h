#ifndef SLAE2_0_FPI_H
#define SLAE2_0_FPI_H

#include <fstream>
#include "../Matrix/CSR.h"
#include "Chebyshev.h"


template<typename T>
std::vector<T> FPI(const CSR<T> &A, const std::vector<T> &b, const std::vector<T> &x0, T tolerance, T tau) {
    std::vector<T> x = x0;
    std::vector<T> x_1 = x;
    std::vector<T> r = A * x - b;
    int count = 1;
    std::ofstream out;
    out.open(R"(C:\Users\user\CLionProjects\SLAE2.0\Data_FPIo.txt)");

    while (!stop_check(A, x_1, b, tolerance)) {
        out << x_1[0] << " " << x_1[3] << "\n";
        x = x_1;
        x_1 = x - (A * x - b) * tau;
        count ++;
    }
    out.close();
    return x_1;
}

template<typename T>
std::vector<T>
FPI_Chebyshev_accelerated(const CSR<T> &A, const std::vector<T> &b, const std::vector<T> &x0, T tolerance, T lam_max,
                          T lam_min, int it) {
    std::vector<T> x = x0;
    std::vector<T> x_1 = x;
    std::vector<T> sol_ = Chebyshev_solutions<T>(it);
    std::vector<int> tau_dist = tau_distrib(it);

    std::vector<T> roots_ = MPI_Cheb_solutions(sol_, tau_dist, lam_max, lam_min);

    int count = 0;

    int i = 1;
    std::vector<double> r;
    std::ofstream out;
    out.open(R"(C:\Users\user\CLionProjects\SLAE2.0\Data_FPIc.txt)");
    out << 0 << " " << 0 << "\n";
    while (!stop_check(A, x_1, b, tolerance)) {
        x = x_1;
        x_1 = x - (A * x - b) * roots_[count];
        count ++;
        if(count == it) count = 0;
        r = A * x - b;
        out << x_1[0] << " " << x_1[3] << "\n";
        i ++;
    }
    out.close();
    return x_1;
}


#endif //SLAE2_0_FPI_H