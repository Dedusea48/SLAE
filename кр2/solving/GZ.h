#ifndef SLAE2_0_GZ_H
#define SLAE2_0_GZ_H
#include <fstream>
#include "../Matrix/CSR.h"


template<typename T>
std::vector<T> Gaus_Zeidel(const CSR<T> &A, const std::vector<T> &b, const std::vector<T> &x0, T tolerance) {
    std::vector<T> x = x0;
    int count = 0;
    while (!stop_check(A, x, b, tolerance)) {

        for (int i = 0; i < x.size(); i++) {
            T sum = b[i];

            for (int j = 0; j < A.get_row(i + 1) - A.get_row(i); j++)
                if (i != j)
                    sum -= A.get_value(A.get_row(i) + j) * x[A.get_col(A.get_row(i) + j)];
            x[i] = sum / A(i, i);
        }
        count++;
    }
    std::cout << count << "\n";
    return x;
}

template<typename T>
std::vector<T> Gaus_Zeidel_Sym(const CSR<T> &A, const std::vector<T> &b, const std::vector<T> &x0, T tolerance) {
    std::vector<T> x = x0;
    while (!stop_check(A, x, b, tolerance)) {
        for (int i = 0; i < x.size(); i++) {
            T sum = b[i];

            for (int j = 0; j < A.get_row(i + 1) - A.get_row(i); j++)
                if (i != j)
                    sum -= A.get_value(A.get_row(i) + j) * x[A.get_col(A.get_row(i) + j)];
            x[i] = sum / A(i, i);
        }

        for (int i = x.size() - 1; i > 0; i--) {
            T sum = b[i];
            for (int j = 0; j < A.get_row(i + 1) - A.get_row(i); j++)
                if (i != j)
                    sum -= A.get_value(A.get_row(i) + j) * x[A.get_col(A.get_row(i) + j)];
            x[i] = sum / A(i, i);
        }
    }
    return x;
}


template<typename T>
std::vector<T> Gaus_Zeidel_Sym_it(const CSR<T> &A, const std::vector<T> &b, const std::vector<T> &x0) {
    std::vector<T> x = x0;
    for (int i = 0; i < x.size(); i++) {
        T sum = b[i];

        for (int j = 0; j < A.get_row(i + 1) - A.get_row(i); j++)
            if (i != j)
                sum -= A.get_value(A.get_row(i) + j) * x[A.get_col(A.get_row(i) + j)];
        x[i] = sum / A(i, i);
    }

    for (int i = x.size() - 1; i >= 0; i--) {
        T sum = b[i];
        for (int j = 0; j < A.get_row(i + 1) - A.get_row(i); j++)
            if (i != j)
                sum -= A.get_value(A.get_row(i) + j) * x[A.get_col(A.get_row(i) + j)];
        x[i] = sum / A(i, i);
    }
    return x;
}

template<typename T>
std::vector<T>
Gaus_Zeidel_accelerated(const CSR<T> &A, const std::vector<T> &b, const std::vector<T> &x, T rho, T tolerance) {
    std::vector<T> x0 = x;
    std::vector<T> x1 = Gaus_Zeidel_Sym_it(A, b, x0);
    std::vector<T> x2 = Gaus_Zeidel_Sym_it(A, b, x1);
    x1 = x2;
    std::vector<double> r = A * x2 - b;

    T mu0 = 1;
    T mu1 = 1 / rho;
    T mu2 = 2 * mu1 / rho - mu0;

    int count = 1;
    std::ofstream out;
    out.open("(C:\Users\user\CLionProjects\SLAE2.0\DResidual(iterations number).txt)");
    out << count << " " << norm(r) << "\n";

    while (!stop_check(A, x2, b, tolerance)) {
        x2 = x1 * (2 * mu1 / (rho * mu2)) - x0 * (mu0 / mu2);
        x0 = x1;
        x1 = Gaus_Zeidel_Sym_it(A, b, x2);
        mu0 = mu1;
        mu1 = mu2;
        mu2 = 2 * mu1 / rho - mu0;
        count++;
        r = A * x2 - b;
        out << count << " " << norm(r) << "\n";
    }
    std::cout << count << "\n";
    return x2;
}


#endif //SLAE2_0_GD_H