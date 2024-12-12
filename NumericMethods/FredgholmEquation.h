#pragma once
#include "stdafx.h"
struct Fredgholm {
    double upper;
    double lower;
    int step;
    double kernel_vlaue;
    std::vector<double> t;
    std::vector<double> nodes;
};

void gauss_legendre(int N, std::vector<double>& nodes, std::vector<double>& weights) {
    nodes.resize(N);
    weights.resize(N);

    if (N == 3) {
        nodes = { -std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0) };
        weights = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
    }
    else if (N == 10) {
        nodes = { -0.9739065285171717, -0.8650633666889845, -0.6794095682990244, -0.4333953941292472, -0.1488743389816312,
                  0.1488743389816312,  0.4333953941292472,  0.6794095682990244,  0.8650633666889845,  0.9739065285171717 };
        weights = { 0.0666713443086881, 0.1494513491505806, 0.2190863625159820, 0.2692667193099963, 0.2955242247147529,
                   0.2955242247147529, 0.2692667193099963, 0.2190863625159820, 0.1494513491505806, 0.0666713443086881 };
    }
}

void transform_to_interval(double a, double b, const std::vector<double>& nodes_ref, const std::vector<double>& weights_ref,
    std::vector<double>& nodes, std::vector<double>& weights) {
    int N = nodes_ref.size();
    nodes.resize(N);
    weights.resize(N);
    double midpoint = 0.5 * (a + b);
    double half_length = 0.5 * (b - a);

    for (int i = 0; i < N; ++i) {
        nodes[i] = midpoint + half_length * nodes_ref[i];
        weights[i] = half_length * weights_ref[i];
    }
}

double kernel_values(double s, double t) {
    return (t / (pow(s, 2)) - 2);
}

double Function(double t) {
    return pow(t, 2) + t / 3 - 1 / 3;
}

std::vector<double> fredgholm_solver_by_trapezium(Fredgholm* equastion_params)
{
    double h = (equastion_params->upper - equastion_params->lower) / equastion_params->step;

    std::vector<double> f_values(equastion_params->step + 1);
    std::vector<double> weights(equastion_params->step + 1);
    equastion_params->t.resize(equastion_params->step + 1);

    std::vector<std::vector<double>> Matrix(
        equastion_params->step + 1,
        std::vector<double>(equastion_params->step + 1)
    );

    std::vector<double> b(equastion_params->step + 1);

    for (size_t i = 0; i <= equastion_params->step; ++i) {
        equastion_params->t[i] = equastion_params->lower + i * h;
        f_values[i] = Function(equastion_params->t[i]);
        weights[i] = (i == 0 || i == equastion_params->step) ? h / 2 : h;
    }

    for (size_t i = 0; i <= equastion_params->step; ++i) {
        for (size_t j = 0; j <= equastion_params->step; ++j) {
            Matrix[i][j] = (i == j ? 1 : 0) - weights[j] * kernel_values(equastion_params->t[j], equastion_params->t[i]);
        }
        b[i] = f_values[i];
    }
    
    return solve_linear_system(Matrix, b);

}

std::vector<double> fregholm_solver_by_gauss(Fredgholm* equastion_params){   

    std::vector<double> nodes, weights;
    std::vector<double> nodes_ref, weights_ref;
    gauss_legendre(equastion_params->step, nodes_ref, weights_ref);
    
    transform_to_interval(equastion_params->lower, equastion_params->upper, nodes_ref, weights_ref, nodes, weights);

    equastion_params->nodes = nodes;

    std::vector<std::vector<double>> A(equastion_params->step, std::vector<double>(equastion_params->step, 0.0));
    std::vector<double> b_vector(equastion_params->step, 0.0);

    for (int i = 0; i < equastion_params->step; ++i) {
        b_vector[i] = Function(nodes[i]);
        for (int j = 0; j < equastion_params->step; ++j) {
            A[i][j] = (i == j ? 1.0 : 0.0) - weights[j] * kernel_values(nodes[j], nodes[i]);
        }
    }

    
    return solveLinearSystem(A, b_vector);
}
