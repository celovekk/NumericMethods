#pragma 
#include <vector>

std::vector<double> solve_linear_system(const std::vector<std::vector<double>>& A, const std::vector<double>& b) {
    int n = A.size();
    std::vector<std::vector<double>> mat = A; 
    std::vector<double> res = b; 
    std::vector<double> x(n);

    // Прямой ход 
    for (int i = 0; i < n; ++i) {
        // Нормализация строки 
        for (int k = i + 1; k < n; ++k) {
            double factor = mat[k][i] / mat[i][i];
            for (int j = i; j < n; ++j)
                mat[k][j] -= factor * mat[i][j];
            res[k] -= factor * res[i];
        }
    }

    // Обратный ход 
    for (int i = n - 1; i >= 0; --i) {
        x[i] = res[i];
        for (int j = i + 1; j < n; ++j)
            x[i] -= mat[i][j] * x[j];
        x[i] /= mat[i][i];
    }
    return x;
}

std::vector<double> solveLinearSystem(std::vector<std::vector<double>>& A, std::vector<double>& b) {
    int N = A.size();
    std::vector<double> solution(N);

    for (int i = 0; i < N; ++i) {
        // Normalize row i
        double diag = A[i][i];
        for (int j = 0; j < N; ++j) {
            A[i][j] /= diag;
        }
        b[i] /= diag;

        // Eliminate column i in subsequent rows
        for (int k = i + 1; k < N; ++k) {
            double factor = A[k][i];
            for (int j = 0; j < N; ++j) {
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
        }
    }

    // Back substitution
    for (int i = N - 1; i >= 0; --i) {
        solution[i] = b[i];
        for (int j = i + 1; j < N; ++j) {
            solution[i] -= A[i][j] * solution[j];
        }
    }

    return solution;
}
