#include <iostream>
#include <cmath>

#include "include/mvn_test.h"

namespace mvn {

mvn_test::mvn_test(Matrix X, const Matrix& S, const Vector& mean) : mahalanobis_centered(X.cols()),
                                                                    mahalanobis_pairwise(X.cols(), X.cols()),
                                                                    fast_impl(mahalanobis_pairwise, 15),
                                                                    wheel(std::random_device()()),
                                                                    the_rest(X.cols()) {
    if (X.cols() == 0 || X.rows() == 0) {
        throw std::invalid_argument("Matrix is empty");
    }
    std::iota(the_rest.begin(), the_rest.end(), 0);
    std::shuffle(the_rest.begin(), the_rest.end(), wheel);
    Eigen::FullPivHouseholderQR<Matrix> qr(S);
    if (!qr.isInvertible()) {
        throw std::logic_error("Non-invertible matrix. Must not happen.");
    }

    Matrix S_inv = qr.inverse();
    n = X.cols();
    p = X.rows();

    if (n <= p) {
        throw std::logic_error("Too few points.");
    }

    for (size_t i = 0; i < n; i++) {
        X.col(i) -= mean;
        const auto& centered = X.col(i);
        mahalanobis_centered[i] = std::exp(-0.25 * centered.transpose() * S_inv * centered);
        for (size_t j = 0; j < n; j++) {
            const auto& lhs = X.col(i);
            const auto& rhs = X.col(j);
            Vector diff = lhs - rhs;
            mahalanobis_pairwise(i, j) = std::exp(-0.5 * diff.transpose() * S_inv * diff);
        }
    }

    for (int i = 0; i < p + 1; i++) {
        add_one();
    }
}

double mvn_test::get_normality_statistic() const {
    if (subsample_size() <= dimensions()) {
        throw std::logic_error("Too few points.");
    }
    double stat = std::pow(3, dimensions() / -2.0);
    stat += (1.0 / ((double) subsample_size() * subsample_size())) * pairwise_stat;
    stat -= (1.0 / subsample_size()) * center_stat * std::pow(2.0, 1 - dimensions() / 2.0);

    return stat;
}

const std::vector<size_t>& mvn_test::current_subset() const {
    return subset;
}

size_t mvn_test::dimensions() const {
    return p;
}

size_t mvn_test::subsample_size() const {
    return subset.size();
}

size_t mvn_test::sample_size() const {
    return n;
}

void mvn_test::swap_once(bool reject_last) {
    if (subset.empty() || the_rest.empty()) {
        throw std::logic_error("Swap points error.");
    }

    if (!reject_last) {
        std::uniform_int_distribution<unsigned> subset_unif(0, subsample_size() - 1);
        std::uniform_int_distribution<unsigned> the_rest_unif(0, the_rest.size() - 1);

        std::swap(subset[subset_unif(wheel)], subset.back());
        std::swap(the_rest[the_rest_unif(wheel)], the_rest.back());
    }

    auto subset_point = subset.back();
    auto replacing_point = the_rest.back();

    fast_impl.unset(subset_point);

    if (!fast_impl.is_effective()) {
        for (auto s: subset) {
            pairwise_stat -= 2 * mahalanobis_pairwise(subset_point, s);
        }

        for (auto s: subset) {
            pairwise_stat += 2 * mahalanobis_pairwise(replacing_point, s);
        }
    } else {
        pairwise_stat -= 2 * fast_impl.score(subset_point);
        pairwise_stat += 2 * fast_impl.score(replacing_point);
    }

    fast_impl.set(replacing_point);

    center_stat -= mahalanobis_centered(subset_point);
    center_stat += mahalanobis_centered(replacing_point);

    std::swap(subset.back(), the_rest.back());
}

void mvn_test::add_one() {
    if (the_rest.empty()) {
        throw std::logic_error("Can't add point to the model.");
    }

    std::uniform_int_distribution<unsigned> random(0, the_rest.size() - 1);
    unsigned pos = random(wheel);
    size_t point = the_rest[pos];
    std::swap(the_rest[pos], the_rest[the_rest.size() - 1]);
    the_rest.pop_back();

    for (size_t s: subset) {
        pairwise_stat += 2 * mahalanobis_pairwise(point, s);
    }

    subset.push_back(point);
    center_stat += mahalanobis_centered(point);
}

}