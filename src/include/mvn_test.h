#ifndef SRC_MVN_TEST_H
#define SRC_MVN_TEST_H

#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

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

enum class mvn_test_method {
    exact,
    rff,
};

mvn_test_method parse_mvn_test_method(const std::string& method);

class mvn_test_base {
public:
    virtual ~mvn_test_base() = default;

    virtual size_t dimensions() const = 0;
    virtual size_t sample_size() const = 0;
    virtual size_t subsample_size() const = 0;

    virtual void add_one() = 0;
    virtual void swap_once(bool reject_last = false) = 0;

    virtual const std::vector<size_t>& current_subset() const = 0;
    virtual double get_normality_statistic() = 0;

    virtual std::shared_ptr<mvn_test_base> clone() const = 0;
};

inline bool operator<(mvn_test_base& lhs, mvn_test_base& rhs) {
    return lhs.get_normality_statistic() < rhs.get_normality_statistic();
}

class mvn_test_exact final : public mvn_test_base {
    class mahalanobis_distances;
    class mvn_stats;

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
    mvn_test_exact(std::shared_ptr<const Matrix> X, const Clustering& clst, const Matrix& cov, const Vector& mean);
    mvn_test_exact(const mvn_test_exact&);

    size_t dimensions() const override;
    size_t sample_size() const override;
    size_t subsample_size() const override;

    void add_one() override;
    void swap_once(bool reject_last = false) override;

    const std::vector<size_t>& current_subset() const override;
    double get_normality_statistic() override;

    std::shared_ptr<mvn_test_base> clone() const override;

private:
    std::vector<double> loglikelihood(const std::vector<int>& ids) const;

    void remove(unsigned i);
    void add(unsigned i);
};

class mvn_test_rff final : public mvn_test_base {
    RandomSampler sampler;

    std::vector<double> pairwise_stat;
    std::vector<double> center_stat;
    std::vector<double> betas;

    struct RFFParams {
        std::shared_ptr<const Matrix> X;
        Matrix A;               // whitening: y = A * x
        Vector mean_whitened;   // A * mean
        Matrix W;               // rff_dim x p
        Vector b;               // rff_dim
        double scale;           // sqrt(2/rff_dim)
    };

    size_t rff_dim;
    std::shared_ptr<const RFFParams> rff;

    std::vector<Eigen::VectorXd> subset_rff_sum; // [beta] : rff_dim
    std::vector<std::unordered_map<size_t, Eigen::VectorXf>> cluster_rff_cache; // [beta]
    std::vector<std::unordered_map<size_t, float>> cluster_center_cache;        // [beta]

    std::shared_ptr<Clustering> clustering;

    size_t p;
    size_t n;

    size_t effect_size;
    int latest_subset_point;

    mutable std::mt19937 wheel;

    std::vector<size_t> subset;

public:
    // rff_dim == 0 means "auto" heuristic.
    mvn_test_rff(std::shared_ptr<const Matrix> X, const Clustering& clst, const Matrix& cov, const Vector& mean, size_t rff_dim);
    mvn_test_rff(const mvn_test_rff&);

    size_t dimensions() const override;
    size_t sample_size() const override;
    size_t subsample_size() const override;

    void add_one() override;
    void swap_once(bool reject_last = false) override;

    const std::vector<size_t>& current_subset() const override;
    double get_normality_statistic() override;

    std::shared_ptr<mvn_test_base> clone() const override;

private:
    void remove(unsigned i);
    void add(unsigned i);
    void ensure_cluster_cached(size_t cluster_id);
};

std::shared_ptr<mvn_test_base> make_mvn_test(std::shared_ptr<const Matrix> X,
                                            const Clustering& clst,
                                            const Vector& mean,
                                            const Matrix& cov,
                                            mvn_test_method method,
                                            size_t rff_dim);

}

#endif // SRC_MVN_TEST_H
