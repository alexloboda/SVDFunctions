#ifndef SRC_MVN_TEST_H
#define SRC_MVN_TEST_H

#include <vector>
#include <memory>
#include <cstdint>
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
    std::shared_ptr<const Matrix> X;

    Matrix S_inv;
    Matrix S_inv_X;
    Vector S_inv_mean;

    std::vector<double> quad_x;
    std::vector<double> x_mu;
    double mu_mu = 0.0;

public:
    mahalanobis_distances(std::shared_ptr<const Matrix> X, const Matrix& cov, const Vector& mean);
    double distance(unsigned i) const;
    double interpoint_distance(unsigned i, unsigned j) const;

    const Matrix& inv_cov() const {
        return S_inv;
    }

    const Matrix& data() const {
        return *X;
    }
};

class mvn_stats {
    std::vector<double> mahalanobis_centered;
    std::vector<std::vector<double>> mahalanobis_pairwise;
    std::vector<std::vector<float>> cluster_features;
    bool feature_mode = false;
    size_t feature_dim = 0;

    const mahalanobis_distances* distances = nullptr;
    const Clustering* clustering = nullptr;
    double k_pw = 0.0;

    double feature_pairwise_stat(size_t i, size_t j) const;
public:
    mvn_stats(const mahalanobis_distances& distances, const Clustering& clst, double beta,
              bool use_nystrom, size_t n_features,
              bool use_rff = false, size_t rff_features = 256, uint32_t seed = 42u);
    mvn_stats() = default;

    bool is_feature_mode() const {
        return feature_mode;
    }

    size_t features_dim() const {
        return feature_dim;
    }

    const float* features_ptr(size_t cluster) const {
        return cluster_features.at(cluster).data();
    }

    double pairwise_stat(size_t i, size_t j) const;
    double sum_pairwise(size_t point, const std::vector<size_t>& ss) const;
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
    std::vector<std::vector<float>> subset_feature_sum;
    std::vector<double> betas;

    std::shared_ptr<Clustering> clustering;

    size_t p;
    size_t n;

    size_t effect_size;
    int latest_subset_point;
    int latest_replacing_point;

    mutable std::mt19937 wheel;

    std::vector<size_t> subset;
    bool use_nystrom = false;
    size_t n_features = 1024;

    bool use_hybrid = false;
    std::shared_ptr<mvn_stats> rff_stats;
    std::vector<size_t> rff_feature_levels;
    std::vector<double> rff_pairwise_stats;
    double rff_center_stat = 0.0;
    std::vector<float> rff_subset_feature_sum;

public:
    mvn_test(std::shared_ptr<const Matrix> X, const Clustering& clst, const Matrix& S, const Vector& mean,
             bool use_nystrom = false, size_t n_features = 1024,
             bool use_hybrid = false, size_t rff_features = 256, uint32_t seed = 42u);
    mvn_test(const mvn_test&);

    size_t dimensions() const;
    size_t sample_size() const;
    size_t subsample_size() const;

    void add_one();
    void swap_once(bool reject_last = false);

    const std::vector<size_t>& current_subset() const;

    double get_normality_statistic();
    double get_normality_statistic() const;

    bool has_aux_statistic() const {
        return use_hybrid && (bool)rff_stats && !rff_feature_levels.empty();
    }

    double get_aux_normality_statistic();
    double get_aux_normality_statistic(size_t level) const;
    size_t aux_statistic_levels() const {
        return rff_feature_levels.size();
    }

    size_t aux_statistic_level_dim(size_t level) const {
        return rff_feature_levels.at(level);
    }

    bool last_swap_has_equal_effect_size() const;
    double exact_delta_last_swap() const;

    friend bool operator<(mvn_test& lhs, mvn_test& rhs);
    std::vector<double> loglikelihood(const std::vector<int>& ids) const;
    std::unique_ptr<mvn_test> clone();

protected:
    void check_aux_state() const;
    void remove(unsigned i);
    void add(unsigned i);

    mvn_test() = default;
};

}

#endif //SRC_MVN_TEST_H
