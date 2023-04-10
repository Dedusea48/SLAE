//
// Created by user on 09.04.2023.
//
#include "iostream"
#include <fstream>
#include <chrono>

#include "Matrix/CSR.h"
#include "solving/FPI.h"
#include "solving/GD.h"
#include "solving/CG.h"
#include "solving/SOR_SSOR.h"


CSR<double> Gen_mat(double a, double b) {
    std::vector<DOK<double>> dok_;
    dok_.reserve(289);
    for (int i = 0; i < 289; i++) {
        dok_.push_back({i, i, 2 * b});
    }

    for (int i = 0; i < 288; i++) {
        dok_.push_back({i + 1, i, a});
        dok_.push_back({i, i + 1, a});
    }

    for (int i = 0; i < 272; i++) {
        dok_.push_back({i + 17, i, a});
        dok_.push_back({i, i + 17, a});
    }

    return CSR<double>{dok_, 289, 289};
}
CSR<double> Gen_mat()
{
    std::vector<DOK<double>> dok_;
    dok_.push_back({0, 0, 3});
    dok_.push_back({1, 1, 3});
    dok_.push_back({2, 2, 4});
    dok_.push_back({3, 3, 4.5});
    return CSR<double>{dok_, 4, 4};
}


int main() {
//Task 1
    double b = 18;
    double a = 6;
    CSR<double> A = Gen_mat(a, b);
    double lam_max = 59.82;
    double lam_min = 12.18;
    std::vector<double> c(289, static_cast<double>(3));
    std::vector<double> x0(289, static_cast<double>(0));
    double t = 1e-13;


    std::vector<double> res_1 = FPI(A, c, x0, t, 1 / lam_max);
    std::vector<double> res_2 = FPI(A, c, x0, t, 2 / (lam_max + lam_min));
    std::vector<double> res_3 = FPI_Chebyshev_accelerated(A, c, x0, t, lam_max, lam_min, 8);
    std::vector<double> res_4 = SSOR(A, c, x0, t, 1.02);

    std::ofstream out;
    out.open(R"(C:\Users\user\PycharmProjects\pythonProject17\speed.txt)");

    double i = -10;
    while (i < 800) {
        auto start = std::chrono::high_resolution_clock::now();

        std::vector<double> res_5 = FPI_Chebyshev_accelerated(A, c, x0, t, lam_max + i, lam_min, 8);

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = duration_cast<std::chrono::microseconds>(stop - start);
        out << duration.count() << "\n";

        i += 0.5;
    }
    out.close();


//Task 2
    std::cout << "Task 2:" << std::endl;

    CSR<double> A2 = Gen_mat();
    std::vector<double> b2(4, 5);
    double c2 = 6;
    std::vector<double> x02(4, static_cast<double>(0));
    double t2 = 1e-13;

    double lam_max2 = 4.5;
    double lam_min2 = 3.0;

    std::vector<double> res2[5];

    res2[0] = FPI(A2, b2, x02, t2, 1.8 / lam_max2);
    res2[1] = FPI(A2, b2, x02, t2, 2 / (lam_max2 + lam_min2));
    res2[2] = Gradient_descent_fast(A2, b2, x02, t2);
    res2[3] = FPI_Chebyshev_accelerated(A2, b2, x02, t2, lam_max2, lam_min2, 8);
    res2[4] = Conjaguate_Gradient(A2, b2, x02, t2);

    double min[5];
    std::cout << "min=";
    for(const auto & item : res2) {
        std::cout << scalar(item, A2 * item) - scalar(item, b2) + c2 << ' ';
    }
}