#include <iostream>
#include <cmath>
#include <unordered_set>

#include "include/mvn_test.h"

namespace mvn {

mvn_test::mvn_test(const Matrix& X, const Clustering& clst, const Matrix& S, const Vector& mean)
        :stats(std::make_shared<mvn_stats>(X, clst, S, mean)),
         pairwise_stat(0.0),
         center_stat(0.0),
         wheel(std::random_device()()),
         the_rest(clst.size()) {
    if (X.cols() == 0 || X.rows() == 0) {
        throw std::invalid_argument("Matrix is empty");
    }
    std::iota(the_rest.begin(), the_rest.end(), 0);
    std::shuffle(the_rest.begin(), the_rest.end(), wheel);

    n = clst.size();
    p = X.rows();
    effect_size = 0;

    if (n <= p) {
        throw std::logic_error("Too few points.");
    }

    while (effect_size < p + 1) {
        add_one();
    }
}

void mvn_stats::compute_full_stats() {
    for (size_t i = 0; i < clustering.size(); i++) {
        double full_row = 0.0;
        for (size_t j = 0; j < clustering.size(); j++) {
            full_row += mahalanobis_pairwise(i, j);
        }
        full_pairwise_stat.push_back(full_row);
    }
}

double mvn_test::get_normality_statistic() const {
    if (effect_size <= dimensions()) {
        throw std::logic_error("Too few points.");
    }
    double stat = std::pow(3, dimensions() / -2.0);
    stat += (1.0 / ((double) effect_size * effect_size)) * pairwise_stat;
    stat -= (1.0 / effect_size) * center_stat * std::pow(2.0, 1 - dimensions() / 2.0);

    return stat;
}

const std::vector<size_t>& mvn_test::current_subset() const {
    return subset;
}

mvn_stats::mvn_stats(const Matrix& X, const Clustering& clst, const Matrix& S, const Vector& mean)
    :mahalanobis_centered(Vector::Zero(clst.size())),
     mahalanobis_pairwise(Matrix::Zero(clst.size(), clst.size())),
     clustering(clst) {
    size_t n = clst.size();

    Eigen::FullPivHouseholderQR<Matrix> qr(S);
    if (!qr.isInvertible()) {
        throw std::logic_error("Non-invertible matrix. Must not happen.");
    }

    Matrix S_inv = qr.inverse();

    for (size_t cl = 0; cl < n; cl++) {
        for (int el: clst.elements(cl)) {
            const auto& centered = X.col(el) - mean;
            mahalanobis_centered[cl] += std::exp(-0.25 * centered.transpose() * S_inv * centered);
            for (size_t pair_cl = 0; pair_cl < n; pair_cl++) {
                for (int pair_el: clst.elements(pair_cl)) {
                    const auto& lhs = X.col(el);
                    const auto& rhs = X.col(pair_el);
                    Vector diff = lhs - rhs;
                    mahalanobis_pairwise(cl, pair_cl) += std::exp(-0.5 * diff.transpose() * S_inv * diff);
                }
            }
        }
    }

    compute_full_stats();
}

size_t mvn_stats::effect_size(size_t i) const {
    return clustering.cluster_size(i);
}

double mvn_stats::pairwise_stat(size_t i, size_t j) const {
    return mahalanobis_pairwise(i, j);
}

double mvn_stats::centered_stat(size_t i) const {
    return mahalanobis_centered(i);
}

double mvn_stats::full_pairwise_statistic(size_t i) const {
    return full_pairwise_stat[i];
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
        throw std::logic_error("Unable to swap points.");
    }

    if (!reject_last) {
        std::uniform_int_distribution<unsigned> subset_unif(0, subsample_size() - 1);
        std::uniform_int_distribution<unsigned> the_rest_unif(0, the_rest.size() - 1);

        std::swap(subset[subset_unif(wheel)], subset.back());
        std::swap(the_rest[the_rest_unif(wheel)], the_rest.back());
    }

    auto subset_point = subset.back();
    auto replacing_point = the_rest.back();

    effect_size = effect_size + stats->effect_size(replacing_point) - stats->effect_size(subset_point);

    if (subset.size() < the_rest.size()) {
        for (auto s: subset) {
            pairwise_stat -= 2 * stats->pairwise_stat(subset_point, s);
        }

        for (auto s: subset) {
            pairwise_stat += 2 * stats->pairwise_stat(replacing_point, s);
        }
    } else {
        pairwise_stat += stats->full_pairwise_statistic(replacing_point);
        pairwise_stat -= stats->full_pairwise_statistic(subset_point);
        for (auto s: the_rest) {
            pairwise_stat += 2 * stats->pairwise_stat(subset_point, s);
        }

        for (auto s: the_rest) {
            pairwise_stat -= 2 * stats->pairwise_stat(replacing_point, s);
        }
    }

    center_stat -= stats->centered_stat(subset_point);
    center_stat += stats->centered_stat(replacing_point);

    std::swap(subset.back(), the_rest.back());
}

void mvn_test::add_one() {
    if (the_rest.empty()) {
        throw std::logic_error("Can't add point to the model.");
    }

    std::uniform_int_distribution<unsigned> random(0, the_rest.size() - 1);
    unsigned pos = random(wheel);
    size_t point = the_rest[pos];
    effect_size += stats->effect_size(point);
    std::swap(the_rest[pos], the_rest[the_rest.size() - 1]);
    the_rest.pop_back();

    for (size_t s: subset) {
        pairwise_stat += 2 * stats->pairwise_stat(point, s);
    }

    subset.push_back(point);
    center_stat += stats->centered_stat(point);
}

void mvn_test::set_solution(std::vector<size_t> solution) {
    std::vector<size_t> whole_set = subset;
    whole_set.insert(whole_set.end(), the_rest.begin(), the_rest.end());
    the_rest.clear();
    subset.clear();
    std::unordered_set<size_t> solution_set(solution.begin(), solution.end());
    for (auto el: whole_set) {
        if (solution_set.find(el) == solution_set.end()) {
            the_rest.push_back(el);
        } else {
            subset.push_back(el);
        }
    }
}

mvn_test::mvn_test(const mvn_test& other)
    :stats(other.stats),
     p(other.p),
     n(other.n),
     effect_size(other.effect_size),
     pairwise_stat(other.pairwise_stat),
     center_stat(other.center_stat),
     wheel{other.wheel()},
     subset(other.subset),
     the_rest(other.the_rest) {}

Clustering::Clustering(const std::vector<int>& clustering) {
    if (clustering.empty()) {
        throw std::invalid_argument("Clustering must not be empty");
    }

    int n_clsuters = *std::max_element(clustering.begin(), clustering.end()) + 1;
    cluster_sizes.resize(n_clsuters);
    clusters.resize(n_clsuters);
    for (size_t i = 0; i < clustering.size(); i++) {
        int cl = clustering[i];
        ++cluster_sizes[cl];
        clusters[cl].push_back(i);
    }
}

size_t Clustering::size() const {
    return clusters.size();
}

std::vector<int> Clustering::elements(size_t i) const {
    return clusters.at(i);
}

size_t Clustering::cluster_size(size_t i) const {
    return clusters.at(i).size();
}

}