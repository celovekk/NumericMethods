#pragma 
#include <vector>

std::vector<double> solve_linear_system(const std::vector<std::vector<double>>& A, const std::vector<double>& b) {
    int n = A.size();
    std::vector<std::vector<double>> mat = A; 
    std::vector<double> res = b; 
    std::vector<double> x(n);

    // ?????? ??? 
    for (int i = 0; i < n; ++i) {
        // ???????????? ?????? 
        for (int k = i + 1; k < n; ++k) {
            double factor = mat[k][i] / mat[i][i];
            for (int j = i; j < n; ++j)
                mat[k][j] -= factor * mat[i][j];
            res[k] -= factor * res[i];
        }
    }

    // ???????? ??? 
    for (int i = n - 1; i >= 0; --i) {
        x[i] = res[i];
        for (int j = i + 1; j < n; ++j)
            x[i] -= mat[i][j] * x[j];
        x[i] /= mat[i][i];
    }
    return x;
}