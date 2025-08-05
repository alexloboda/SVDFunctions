#ifndef SRC_KRON_H
#define SRC_KRON_H

#include "third-party/cxxpool.h"
#include "third-party/zstr/zstr.hpp"

#include <vector>
#include <fstream>
#include <memory>

#include <RcppEigen.h>

namespace matching {

namespace impl {

using matrix_t = Eigen::MatrixXd;

class one_spot_approximation {
    int k;
    int m;
    std::vector<matrix_t> matrices;
public:
    one_spot_approximation(const matrix_t& matrix, int m, int k);
    one_spot_approximation(int m, int k);
    matrix_t get_approximation() const;
    double calculate(const matrix_t& sigma, double c_e) const;
    friend std::ostream& operator<<(std::ostream& os, const one_spot_approximation& spot);
    friend std::istream& operator>>(std::istream& is, one_spot_approximation& spot);
};

class one_degree_approximation {
    int m;
    int k;
    std::vector<one_spot_approximation> spots;
public:
    one_degree_approximation(const matrix_t& outer, int m, int k);
    one_degree_approximation(int m, int k);
    double calculate(const matrix_t& sigma, double c_e) const;
    int n_spots() const;
    friend std::ostream& operator<<(std::ostream& os, const one_degree_approximation& spot);
    friend std::istream& operator>>(std::istream& is, one_degree_approximation& spot);
};

class kronecker_approximation {
    std::vector<one_degree_approximation> degrees;
    int m;
public:
    kronecker_approximation(const matrix_t& A, const matrix_t& B, int max_degree);
    kronecker_approximation();
    double compression() const;
    double calculate(const matrix_t& sigma, double c_e) const;
    friend std::ostream& operator<<(std::ostream& os, const kronecker_approximation& spot);
    friend std::istream& operator>>(std::istream& is, kronecker_approximation& spot);
};


}

class kronecker_calculator {
    zstr::ifstream fin;
public:
    kronecker_calculator(std::string filename);
    std::vector<std::vector<double>> calculate(Eigen::MatrixXd sigma, double c_e);
};

class kronecker_preprocessor {
    std::shared_ptr<Eigen::MatrixXd> matrix;
    std::vector<std::vector<int>> clusters;

    zstr::ofstream fout;
    std::vector<std::future<std::unique_ptr<impl::kronecker_approximation>>> futures;

public: 
    kronecker_preprocessor(std::shared_ptr<Eigen::MatrixXd> matrix, 
                           std::vector<std::vector<int>> clusters, 
                           std::string filename);

    void process(unsigned threads, unsigned batch_size, unsigned max_degree);
    double write_futures();
private:

};
}

#endif 