#pragma once
#include "stdafx.h"
struct Fredgholm {
    double upper;
    double lower;
    int step;
    double kernel_vlaue;
};

double kernel_values(double s, double t) {
    return (t / (pow(s, 2)) - 2);
}

double Function(double t) {
    return pow(t, 2) + t / 3 - 1 / 3;
}

std::vector<double> fredgholm_solver(Fredgholm* equastion_params)
{
    double h = (equastion_params->upper - equastion_params->lower) / equastion_params->step;

    std::vector<double> t(equastion_params->step + 1);
    std::vector<double> weights(equastion_params->step + 1);
    std::vector<double> f_values(equastion_params->step + 1);

    std::vector<std::vector<double>> Matrix(
        equastion_params->step + 1,
        std::vector<double>(equastion_params->step + 1)
    );

    std::vector<double> b(equastion_params->step + 1);

    for (size_t i = 0; i <= equastion_params->step; ++i) {
        t[i] = equastion_params->lower + i * h;
        f_values[i] = Function(t[i]);
        weights[i] = (i == 0 || i == equastion_params->step) ? h / 2 : h;
    }

    for (size_t i = 0; i <= equastion_params->step; ++i) {
        for (size_t j = 0; j <= equastion_params->step; ++j) {
            Matrix[i][j] = (i == j ? 1 : 0) - weights[j] * kernel_values(t[j], t[i]);
        }
        b[i] = f_values[i];
    }
    return solve_linear_system(Matrix, b);

}