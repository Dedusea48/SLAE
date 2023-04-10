#ifndef SLAE2_0_GD_H
#define SLAE2_0_GD_H

#include <fstream>
#include "../Matrix/CSR.h"


template<typename T>
std::vector<T> Gradient_descent_fast(const CSR<T> &A, const std::vector<T> &b, const std::vector<T> &x0, T tolerance)
{
    std::vector<T> x = x0;
    std::vector<T> r = A * x - b;
    T a;
    std::ofstream out;
    out.open(R"(C:\Users\user\CLionProjects\SLAE2.0\Data_GD.txt)");
    out << 0 << " " << 0 << "\n";
    while(!stop_check(A, x, b, tolerance))
    {
        a = scalar(r, r) / scalar(r, A * r);
        x = x - r * a;
        r = A * x - b;
        out << x[0] << " " << x[3] << "\n";
    }
    out.close();
    return x;
}

#endif //SLAE2_0_GD_H
