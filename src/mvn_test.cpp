#include <iostream>
#include <cmath>
#include <unordered_set>

#include "include/mvn_test.h"

namespace mvn {

mvn_test::mvn_test(const Matrix& X, const Clustering& clst, const Matrix& S, const Vector& mean)
        :clustering(std::make_shared<Clustering>(clst)),
         wheel(std::random_device()()),
         the_rest(clst.size()) {
    if (X.cols() == 0 || X.rows() == 0) {
        throw std::invalid_argument("Matrix is empty");
    }
    double beta = 0.2;
    while (beta < 2) {
        stats.push_back(std::make_shared<mvn_stats>(X, clst, S, mean, beta));
        beta *= 2;
    }
    pairwise_stat.resize(stats.size(), 0.0);
    center_stat.resize(stats.size(), 0.0);

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

double mvn_test::get_normality_statistic() const {
    if (effect_size <= dimensions()) {
        throw std::logic_error("Too few points.");
    }
    double max_stat = 0.0;
    for (size_t i = 0; i < stats.size(); i++) {
        double stat = std::pow(1 + 2 * std::pow(stats[i]->beta(), 2), dimensions() / -2.0);
        stat += (1.0 / ((double) effect_size * effect_size)) * pairwise_stat[i];
        stat -= (2.0 / (effect_size * std::pow(1 + std::pow(stats[i]->beta(), 2.0), dimensions() / 2.0))) * center_stat[i];
        max_stat = std::max(max_stat, stat);
    }

    return max_stat;
}

const std::vector<size_t>& mvn_test::current_subset() const {
    return subset;
}

mvn_stats::mvn_stats(const Matrix& X, const Clustering& clst, const Matrix& S, const Vector& mean, double beta)
    :mahalanobis_centered(Vector::Zero(clst.size())),
     mahalanobis_pairwise(Matrix::Zero(clst.size(), clst.size())),
     beta_val(beta) {
    size_t n = clst.size();

    Eigen::FullPivHouseholderQR<Matrix> qr(S);
    if (!qr.isInvertible()) {
        throw std::logic_error("Non-invertible matrix. Must not happen.");
    }

    Matrix S_inv = qr.inverse();

    double k_center = -(beta * beta / (2 * (1 + beta * beta)));
    double k_pw = -(beta * beta) / 2.0;

    Matrix ximu = X.transpose() * S_inv * mean;
    Matrix muxi = mean.transpose() * S_inv * X;
    Matrix xixj = X.transpose() * S_inv * X;
    double mumu = mean.transpose() * S_inv * mean;

    for (size_t cl = 0; cl < n; cl++) {
        for (int el: clst.elements(cl)) {
            double distance = xixj(el, el) - ximu(el) - muxi(el) + mumu;
            mahalanobis_centered[cl] += std::exp(k_center * distance);
            for (size_t pair_cl = 0; pair_cl < n; pair_cl++) {
                for (int pair_el: clst.elements(pair_cl)) {
                    double pair_dist = xixj(el, el) - xixj(el, pair_el) - xixj(pair_el, el) + xixj(pair_el, pair_el);
                    if (cl == pair_cl) {
                        mahalanobis_pairwise(cl, pair_cl) += 0.5 * std::exp(k_pw * pair_dist);
                    } else {
                        mahalanobis_pairwise(cl, pair_cl) += std::exp(k_pw * pair_dist);
                    }
                }
            }
        }
    }
}

double mvn_stats::pairwise_stat(size_t i, size_t j) const {
    return mahalanobis_pairwise(i, j);
}

double mvn_stats::centered_stat(size_t i) const {
    return mahalanobis_centered(i);
}

double mvn_stats::beta() const {
    return beta_val;
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


    effect_size += clustering->cluster_size(replacing_point) - clustering->cluster_size(subset_point);

    for (size_t i = 0; i < stats.size(); i++) {
        for (auto s: subset) {
            pairwise_stat[i] -= 2 * stats[i]->pairwise_stat(subset_point, s);
        }
    }

    std::swap(subset.back(), the_rest.back());

    for (size_t i = 0; i < stats.size(); i++) {
        for (auto s: subset) {
            pairwise_stat[i] += 2 * stats[i]->pairwise_stat(replacing_point, s);
        }
    }

    for (size_t i = 0; i < stats.size(); i++) {
        center_stat[i] -= stats[i]->centered_stat(subset_point);
        center_stat[i] += stats[i]->centered_stat(replacing_point);
    }
}

void mvn_test::add_one() {
    if (the_rest.empty()) {
        throw std::logic_error("Can't add point to the model.");
    }

    std::uniform_int_distribution<unsigned> random(0, the_rest.size() - 1);
    unsigned pos = random(wheel);
    size_t point = the_rest[pos];
    effect_size += clustering->cluster_size(point);
    std::swap(the_rest[pos], the_rest[the_rest.size() - 1]);
    the_rest.pop_back();

    subset.push_back(point);
    for (size_t i = 0; i < stats.size(); i++) {
        for (size_t s: subset) {
            pairwise_stat[i] += 2 * stats[i]->pairwise_stat(point, s);
        }

        center_stat[i] += stats[i]->centered_stat(point);
    }
}

mvn_test& mvn_test::operator=(const mvn_test& other) {
    mvn_test copy(other);
    swap(copy);
    return *this;
}

bool operator<(const mvn_test& lhs, const mvn_test& rhs) {
    return lhs.get_normality_statistic() < rhs.get_normality_statistic();
}

mvn_test::mvn_test(const mvn_test& other)
    :stats(other.stats),
     clustering(other.clustering),
     p(other.p),
     n(other.n),
     effect_size(other.effect_size),
     pairwise_stat(other.pairwise_stat),
     center_stat(other.center_stat),
     wheel{other.wheel()},
     subset(other.subset),
     the_rest(other.the_rest) {}

void mvn_test::swap(mvn_test& other) {
    std::swap(other.stats, stats);
    std::swap(other.clustering, clustering);
    std::swap(other.p, p);
    std::swap(other.n, n);
    std::swap(other.subset, subset);
    std::swap(other.the_rest, the_rest);
    std::swap(other.wheel, wheel);
    std::swap(other.pairwise_stat, pairwise_stat);
    std::swap(other.effect_size, effect_size);
    std::swap(other.center_stat, center_stat);
}

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