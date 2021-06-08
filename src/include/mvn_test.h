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
    std::vector<double> pairwise_stat;
    std::vector<double> center_stat;
    std::vector<double> betas;

    std::shared_ptr<Clustering> clustering;

    size_t p;
    size_t n;
    size_t effect_size;

    mutable std::mt19937 wheel;

    std::vector<size_t> subset;
    std::vector<size_t> the_rest;

public:
    size_t dimensions() const;
    size_t sample_size() const;
    size_t subsample_size() const;

    void add_one();
    void swap_once(bool reject_last = false);

    const std::vector<size_t>& current_subset() const;

    double get_normality_statistic();

    friend bool operator<(mvn_test& lhs, mvn_test& rhs);
    virtual std::vector<double> loglikelihood(const std::vector<int>& ids) const = 0;
    virtual std::unique_ptr<mvn_test> clone() = 0;

protected:
    virtual void remove(unsigned i) = 0;
    virtual void add(unsigned i) = 0;
    virtual void compute_statistics() = 0;

    mvn_test(const mvn_test&);
    mvn_test(std::shared_ptr<const Matrix> X, const Clustering& clst);
    mvn_test() = default;
};

class mvn_test_fixed: public mvn_test {
    std::shared_ptr<mahalanobis_distances> distances;
    std::vector<std::shared_ptr<mvn_stats>> stats;
public:
    mvn_test_fixed(std::shared_ptr<const Matrix> X, const Clustering& clst, const Matrix& S, const Vector& mean);
    mvn_test_fixed(const mvn_test_fixed& other);

    void compute_statistics() override;
    void remove(unsigned i) override;
    void add(unsigned i) override;
    std::unique_ptr<mvn_test> clone() override;
    std::vector<double> loglikelihood(const std::vector<int>& ids) const override;
};

class mvn_test_gen: public mvn_test {
    std::shared_ptr<const Matrix> X;
public:
    mvn_test_gen(std::shared_ptr<const Matrix> X, const Clustering& clst);
    mvn_test_gen(const mvn_test_gen& other);

    void compute_statistics() override;
    void remove(unsigned i) override;
    void add(unsigned i) override;
    std::unique_ptr<mvn_test> clone() override;
    std::vector<double> loglikelihood(const std::vector<int>& ids) const override;
};

}

#endif //SRC_MVN_TEST_H
