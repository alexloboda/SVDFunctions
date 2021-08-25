#ifndef SRC_MVN_TEST_H
#define SRC_MVN_TEST_H

#include <vector>
#include <memory>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <random>

namespace mvn {

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

class RandomSampler {
    std::uniform_real_distribution<double> runif;
    mutable std::mt19937 wheel;
    std::vector<double> original;
    std::vector<double> segment_tree;
    std::vector<size_t> active_tree;

    size_t size;
public:
    RandomSampler();
    RandomSampler(const std::vector<double>& logscale, long seed);
    RandomSampler(RandomSampler&&) = default;
    RandomSampler(const RandomSampler& other);
    RandomSampler& operator=(RandomSampler&&);

    void disable(size_t n);
    void enable(size_t n);
    size_t sample();
    size_t n_active() const;
private:
    std::pair<size_t, size_t> children(size_t node) const;
    static bool is_root(size_t node);
    bool is_leaf(size_t node) const;
    bool is_active(size_t node) const;
    size_t el_pos(size_t el) const;
    static size_t parent(size_t node);

    static double sum_log(double l, double r);

    void update_inner_node(size_t node);
    void update(size_t node);
};

class Clustering {
    std::vector<int> cluster_sizes;
    std::vector<std::vector<int>> clusters;
public:
    explicit Clustering(const std::vector<int>& clustering);
    Clustering() = default;
    const std::vector<int>& elements(size_t i) const;
    size_t size() const;
    size_t cluster_size(size_t i) const;
};

class mahalanobis_distances {
    Matrix distances;
    Vector diag;
    Matrix muxi;
    Matrix ximu;
    double mumu;
public:
    mahalanobis_distances(std::shared_ptr<const Matrix> X, const Matrix& cov, const Vector& mean);
    double distance(unsigned i) const;
    double interpoint_distance(unsigned i, unsigned j) const;
};

class mvn_stats {
    Vector mahalanobis_centered;
    Matrix mahalanobis_pairwise;
public:
    mvn_stats(const mahalanobis_distances& distances, const Clustering& clst, double beta);
    mvn_stats() = default;

    double pairwise_stat(size_t i, size_t j) const;
    double centered_stat(size_t i) const;
private:
};

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
    void remove(unsigned i);
    void add(unsigned i);

    mvn_test() = default;
};

}

#endif //SRC_MVN_TEST_H
