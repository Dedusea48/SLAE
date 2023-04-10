//
// Created by user on 09.04.2023.
//

#ifndef SLAE2_0_CHEBYSHEV_H
#define SLAE2_0_CHEBYSHEV_H


#include "../Matrix/CSR.h"

template<typename T>
std::vector<T> Chebyshev_solutions(int n) {
    std::vector<T> solutions_(n);
    T sina = std::sin(M_PI / static_cast<T>(n));
    T cosa = std::cos(M_PI / static_cast<T>(n));
    T sinb = std::sin(M_PI / (2 * static_cast<T>(n)));

    solutions_[0] = std::cos(M_PI / (2 * static_cast<T>(n)));
    for (int i = 1; i < n; i++) {
        solutions_[i] = solutions_[i - 1] * cosa - sinb * sina;
        sinb = solutions_[i - 1] * sina + sinb * cosa;
    }
    //std::cout << solutions_;
    return solutions_;
}

std::vector<int> tau_distrib(int n) {
    std::vector<int> order_(n);
    order_[0] = 1;
    for (int i = 1; i < n; i *= 2) {
        for (int j = 0; j < i; j+=1) {
            order_[n * (j * 2 + 1) / (2 * i)] = 2 * i  + 1 - order_[n * 2 * j / (2 * i)];
        }
    }
    return order_;
}

template<typename T>
std::vector<T> MPI_Cheb_solutions(const std::vector<T>& sol_, const std::vector<int>& order_, T lam_max, T lam_min) {
    std::vector<T> res_(sol_.size());
    for(int i = 0; i < sol_.size(); i ++)
    {
        res_[i] = 2 / ((lam_max + lam_min) + (lam_max - lam_min) * sol_[order_[i] - 1]);
    }
    return res_;
}



#endif //SLAE2_0_CHEBYSHEV_H
