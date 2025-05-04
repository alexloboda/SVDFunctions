#ifndef SRC_MVN_TEST_H
#define SRC_MVN_TEST_H

#include <vector>
#include <memory>
#include <optional>
#include <RcppEigen.h>
#include "third-party/zstr/zstr.hpp"
#include "mvn_stats.h"
#include "mvn_clst.h" 
#include "RandomSampler.h"

// [[Rcpp::depends(RcppEigen)]]

#include <random>

namespace mvn {

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

class mvn_test {
protected:
    std::shared_ptr<mahalanobis_distances> distances;
    std::vector<std::shared_ptr<mvn_stats>> stats;
    RandomSampler sampler;

    std::vector<double> pairwise_stat;
    std::vector<double> center_stat;
    std::vector<double> betas;

    std::shared_ptr<Clustering> clustering;

    size_t p;
    size_t n;

    size_t effect_size;
    int latest_subset_point;

    mutable std::mt19937 wheel;

    std::vector<size_t> subset;

public:
    mvn_test(std::shared_ptr<const Matrix> X, const Clustering& clst, const Matrix& S, const Vector& mean);
    mvn_test(std::shared_ptr<const Matrix> X, const Clustering& clst, const Matrix& S, const Vector& mean, const std::string& filename);
    mvn_test(const mvn_test&);

    size_t dimensions() const;
    size_t sample_size() const;
    size_t subsample_size() const;

    void add_one();
    void swap_once(bool reject_last = false);

    const std::vector<size_t>& current_subset() const;

    double get_normality_statistic();

    friend bool operator<(mvn_test& lhs, mvn_test& rhs);
    std::vector<double> loglikelihood(const std::vector<int>& ids) const;
    std::unique_ptr<mvn_test> clone();

protected:
    void remove(unsigned point); // Updated parameter name for clarity
    void add(unsigned point);    // Updated parameter name for clarity

    mvn_test() = default;

private:
    void initialize_common(std::shared_ptr<const Matrix> X, const Clustering& clst, const Matrix& S, const Vector& mean, std::optional<std::string> filename);
    void initialize_stats(const Clustering& clst, std::optional<std::string> filename); // Updated method name for consistency
    std::vector<double> compute_loglikelihoods(const Clustering& clst) const;          // Added missing declaration
    void update_stats(unsigned point, bool is_addition);                               // Added missing declaration
    void validate_input(std::shared_ptr<const Matrix> X, const Clustering& clst);
    void initialize_sampler(const Clustering& clst);
};

}

#endif //SRC_MVN_TEST_H
