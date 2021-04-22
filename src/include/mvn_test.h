#ifndef SRC_MVN_TEST_H
#define SRC_MVN_TEST_H

#include <vector>
#include <memory>
#include <eigen3/Eigen/Dense>
#include <random>

namespace mvn {

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

class Clustering {
    std::vector<int> cluster_sizes;
    std::vector<std::vector<int>> clusters;
public:
    explicit Clustering(const std::vector<int>& clustering);
    Clustering() = default;
    std::vector<int> elements(size_t i) const;
    size_t size() const;
    size_t cluster_size(size_t i) const;
};

class mvn_stats {
    Vector mahalanobis_centered;
    Matrix mahalanobis_pairwise;

    double beta_val;
public:
    mvn_stats(const Matrix& X, const Clustering& clst, const Matrix& S, const Vector& mean, double beta);
    mvn_stats() = default;

    double pairwise_stat(size_t i, size_t j) const;
    double centered_stat(size_t i) const;
    double beta() const;
private:
};

class mvn_test {
    std::vector<std::shared_ptr<mvn_stats>> stats;
    std::shared_ptr<Clustering> clustering;

    size_t p;
    size_t n;
    size_t effect_size;

    std::vector<double> pairwise_stat;
    std::vector<double> center_stat;

    mutable std::mt19937 wheel;

    std::vector<size_t> subset;
    std::vector<size_t> the_rest;

public:
    mvn_test() = default;
    mvn_test(const mvn_test&);
    mvn_test& operator=(const mvn_test&);
    mvn_test(const Matrix& X, const Clustering& clst, const Matrix& S, const Vector& mean);
    size_t dimensions() const;
    size_t sample_size() const;
    size_t subsample_size() const;

    void add_one();
    void swap_once(bool reject_last = false);
    void swap(mvn_test& other);

    const std::vector<size_t>& current_subset() const;
    double get_normality_statistic() const;

    double get_pairwise_statistic() const;
    double get_centered_statistic() const;

    friend bool operator<(const mvn_test& lhs, const mvn_test& rhs);
};

}

#endif //SRC_MVN_TEST_H
