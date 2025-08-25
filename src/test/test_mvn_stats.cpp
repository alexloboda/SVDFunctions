// Minimal smoke test harness for mvn_stats + kronecker algo build
#include "../include/kronecker.h"
#include <Eigen/Dense>
#include <iostream>

int main() {
    try {
        Eigen::MatrixXd A(4, 3); A.setRandom();
        Eigen::MatrixXd B(5, 3); B.setRandom();
        matching::impl::tpm_params params;
        params.sample_cap = 100; // small sample
        params.max_components = 3;
        params.max_iters = 20;
        params.tol = 1e-5;
        params.seed = 42;
        params.use_irls = true;
        params.enforce_identity_scalar = true;

        matching::impl::kronecker_approximation approx(A, B, 2, params);
        Eigen::MatrixXd S = Eigen::MatrixXd::Identity(3, 3);
        double val = approx.calculate(S, 0.5);
        std::cout << "score=" << val << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
