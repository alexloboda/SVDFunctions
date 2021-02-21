#ifndef SRC_MVN_TEST_H
#define SRC_MVN_TEST_H

#include <vector>
#include <eigen3/Eigen/Dense>
#include <random>

namespace mvn {

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

class batch {
    std::vector<double> scores;
public:
    batch(const mvn::Vector& v, unsigned start, unsigned stop);
    double get_batch_score (unsigned state) const;
};

class state {
    std::vector<unsigned> bits;
    size_t chunk_size;
    unsigned subset_size;
public:
    state() = default;
    state(const state&) = default;
    state(size_t size, size_t chunk_size);

    void set(unsigned i) {
        subset_size++;
        bits.at(i / chunk_size) = bits.at(i / chunk_size) | (1 << i % chunk_size);
    }

    void unset(unsigned i) {
        subset_size--;
        bits.at(i / chunk_size) = bits.at(i / chunk_size) & (~(1 << i % chunk_size));
    }

    unsigned effective_size() const {
        return subset_size;
    }

    std::vector<unsigned> get_states() const {
        return bits;
    }
};

class four_russians {
    size_t n;
    std::vector<std::vector<batch>> batches;
    state current;
public:
    four_russians() = default;
    four_russians(const four_russians&) = default;
    four_russians(const Matrix& pairwise_mahalanobis, size_t batch_size);

    void set(unsigned i);
    void unset(unsigned i);

    double score(unsigned row) const;
    bool is_effective() const;
};


class mvn_test {
    Vector mahalanobis_centered;
    Matrix mahalanobis_pairwise;

    four_russians fast_impl;

    size_t p;
    size_t n;

    double pairwise_stat;
    double center_stat;

    std::mt19937 wheel;

    std::vector<size_t> subset;
    std::vector<size_t> the_rest;

public:
    mvn_test() = default;
    mvn_test(const mvn_test&) = default;
    mvn_test(Matrix X, const Matrix& S, const Vector& mean);
    size_t dimensions() const;
    size_t sample_size() const;
    size_t subsample_size() const;

    void add_one();
    void swap_once(bool reject_last = false);

    const std::vector<size_t>& current_subset() const;
    double get_normality_statistic() const;
};

}

#endif //SRC_MVN_TEST_H
