//
// Created by user on 09.04.2023.
//

#ifndef SLAE2_0_SOR_SSOR_H
#define SLAE2_0_SOR_SSOR_H


#include <fstream>
#include "../Matrix/CSR.h"

template<typename T>
std::vector<T> SOR(const CSR<T> &A, const std::vector<T> &b, const std::vector<T> &x0, T tolerance, T w) {
    std::vector<T> x = x0;
    while (!stop_check(A, x, b, tolerance)) {
        for (int i = 0; i < x.size(); i++) {
            T sum = b[i];
            for (int j = 0; j < A.get_row(i + 1) - A.get_row(i); j++)
                if (i != A.get_col(A.get_row(i) + j))
                    sum -= A.get_value(A.get_row(i) + j) * x[A.get_col(A.get_row(i) + j)];
            x[i] = (1 - w) * x[i] + w * sum / A(i, i);
        }
    }
    return x;
}


template<typename T>
std::vector<T> SSOR(const CSR<T> &A, const std::vector<T> &b, const std::vector<T> &x0, T tolerance, T w) {
    std::vector<T> x = x0;

    std::ofstream out;
    out.open(R"(C:\Users\user\CLionProjects\SLAE2.0\Data_SOR_SSOR.txt)");
    std::vector<T> r = A * x - b;
    int count = 1;

    while (!stop_check(A, x, b, tolerance)) {
        out << count << " " << norm(r) << "\n";
        for (int i = 0; i < x.size(); i++) {
            T sum = b[i];
            for (int j = 0; j < A.get_row(i + 1) - A.get_row(i); j++)
                if (i != A.get_col(A.get_row(i) + j))
                    sum -= A.get_value(A.get_row(i) + j) * x[A.get_col(A.get_row(i) + j)];
            x[i] = (1 - w) * x[i] + w * sum / A(i, i);
        }

        for (int i = x.size() - 1; i > 0; i--) {
            T sum = b[i];
            for (int j = 0; j < A.get_row(i + 1) - A.get_row(i); j++)
                if (i != A.get_col(A.get_row(i) + j))
                    sum -= A.get_value(A.get_row(i) + j) * x[A.get_col(A.get_row(i) + j)];
            x[i] = (1 - w) * x[i] + w * sum / A(i, i);
        }
        count ++;
        r = A * x - b;
    }
    out.close();
    return x;
}

#endif //SLAE2_0_SOR_SSOR_H
