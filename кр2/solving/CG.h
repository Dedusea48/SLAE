//
// Created by user on 09.04.2023.
//

#ifndef SLAE2_0_CG_H
#define SLAE2_0_CG_H


#include <fstream>
#include "../Matrix/CSR.h"

template<typename T>
std::vector<T> Conjaguate_Gradient(const CSR<T> &A, const std::vector<T> &b, const std::vector<T> &x0, T tolerance)
{
    std::vector<T> x = x0;
    std::vector<T> r = A * x - b, r1;
    std::vector<T> d = r;
    T a, beta;
    int count = 0;

    std::ofstream out;
    out.open(R"(C:\Users\user\CLionProjects\SLAE2.0\Data_CG.txt)");
    out << 0 << " " << 0 << "\n";

    while(!stop_check(A, x, b, tolerance))
    {
        a = scalar(d, r) / scalar(d, A * d);
        x = x - d * a;
        r1 = A * x - b;
        beta = scalar(r1, r1) / scalar(r, r);
        r = r1;
        d = r1 + d * beta;
        count ++;
        out << x[0] << " " << x[3] << "\n";
    }
    return x;
}


#endif //SLAE2_0_CG_H
