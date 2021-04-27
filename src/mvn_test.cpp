#include <iostream>
#include <cmath>
#include <unordered_set>

#include "include/mvn_test.h"

namespace mvn {

mvn_test::mvn_test(std::shared_ptr<const Matrix> X, const Clustering& clst)
        :clustering(std::make_shared<Clustering>(clst)),
         wheel(std::random_device()()),
         the_rest(clst.size()) {
    if (X->cols() == 0 || X->rows() == 0) {
        throw std::invalid_argument("Matrix is empty");
    }

    double beta = 0.2;
    while (beta < 2) {
        betas.push_back(beta);
        beta *= 2;
    }
    pairwise_stat.resize(betas.size(), 0.0);
    center_stat.resize(betas.size(), 0.0);

    std::iota(the_rest.begin(), the_rest.end(), 0);
    std::shuffle(the_rest.begin(), the_rest.end(), wheel);

    n = clst.size();
    p = X->rows();
    effect_size = 0;

    if (n <= p) {
        throw std::logic_error("Too few points.");
    }
}

double mvn_test::get_normality_statistic() const {
    if (effect_size <= dimensions()) {
        throw std::logic_error("Too few points.");
    }
    double max_stat = 0.0;
    for (size_t i = 0; i < betas.size(); i++) {
        double stat = std::pow(1 + 2 * std::pow(betas[i], 2), dimensions() / -2.0);
        stat += (1.0 / ((double) effect_size * effect_size)) * pairwise_stat[i];
        stat -= (2.0 / (effect_size * std::pow(1 + std::pow(betas[i], 2.0), dimensions() / 2.0))) * center_stat[i];
        max_stat = std::max(max_stat, stat);
    }

    return max_stat;
}

const std::vector<size_t>& mvn_test::current_subset() const {
    return subset;
}

mvn_stats::mvn_stats(const mahalanobis_distances& distances, const Clustering& clst, double beta)
    :mahalanobis_centered(Vector::Zero(clst.size())),
     mahalanobis_pairwise(Matrix::Zero(clst.size(), clst.size())),
     beta_val(beta) {
    size_t n = clst.size();

    double k_center = -(beta * beta / (2 * (1 + beta * beta)));
    double k_pw = -(beta * beta) / 2.0;

    for (size_t cl = 0; cl < n; cl++) {
        for (int el: clst.elements(cl)) {
            mahalanobis_centered[cl] += std::exp(k_center * distances.distance(el));
            for (size_t pair_cl = 0; pair_cl < n; pair_cl++) {
                for (int pair_el: clst.elements(pair_cl)) {
                    if (cl == pair_cl) {
                        mahalanobis_pairwise(cl, pair_cl) += 0.5 * std::exp(k_pw * distances.interpoint_distance(el, pair_el));
                    } else {
                        mahalanobis_pairwise(cl, pair_cl) += std::exp(k_pw * distances.interpoint_distance(el, pair_el));
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

    remove(subset_point);

    std::swap(subset.back(), the_rest.back());

    add(replacing_point);

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
    add(point);
}

bool operator<(const mvn_test& lhs, const mvn_test& rhs) {
    return lhs.get_normality_statistic() < rhs.get_normality_statistic();
}

mvn_test::mvn_test(const mvn_test& other)
    :clustering(other.clustering),
     p(other.p),
     n(other.n),
     effect_size(other.effect_size),
     wheel{other.wheel()},
     subset(other.subset),
     the_rest(other.the_rest) {}

void mvn_test::swap(mvn_test& other) {
    std::swap(other.clustering, clustering);
    std::swap(other.p, p);
    std::swap(other.n, n);
    std::swap(other.subset, subset);
    std::swap(other.the_rest, the_rest);
    std::swap(other.wheel, wheel);
    std::swap(other.effect_size, effect_size);
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

mahalanobis_distances::mahalanobis_distances(std::shared_ptr<const Matrix> X, const Matrix& S, const Vector& mean) {
    Eigen::FullPivHouseholderQR<Matrix> qr(S);
    if (!qr.isInvertible()) {
        throw std::logic_error("Non-invertible matrix. Must not happen.");
    }

    Matrix S_inv = qr.inverse();

    ximu = X->transpose() * S_inv * mean;
    muxi = mean.transpose() * S_inv * *X;
    xixj = X->transpose() * S_inv * *X;
    mumu = mean.transpose() * S_inv * mean;
}

double mahalanobis_distances::interpoint_distance(unsigned i, unsigned j) const {
    return xixj(i, i) - xixj(i, j) - xixj(j, i) + xixj(j, j);
}

double mahalanobis_distances::distance(unsigned el) const {
    return xixj(el, el) - ximu(el) - muxi(el) + mumu;
}

mvn_test_fixed::mvn_test_fixed(std::shared_ptr<const Matrix> X, const Clustering& clst, const Matrix& S, const Vector& mean)
        :mvn_test(X, clst) {
    mahalanobis_distances distances(X, S, mean);

    for (double beta: betas) {
        stats.push_back(std::make_shared<mvn_stats>(distances, clst, beta));
    }

    while (effect_size < p + 1) {
        add_one();
    }
}

void mvn_test_fixed::remove(unsigned point) {
    for (size_t i = 0; i < stats.size(); i++) {
        for (auto s: subset) {
            pairwise_stat[i] -= 2 * stats[i]->pairwise_stat(point, s);
        }
    }

    for (size_t i = 0; i < stats.size(); i++) {
        center_stat[i] -= stats[i]->centered_stat(point);
    }
}

void mvn_test_fixed::add(unsigned int point) {
    for (size_t i = 0; i < stats.size(); i++) {
        for (auto s: subset) {
            pairwise_stat[i] += 2 * stats[i]->pairwise_stat(point, s);
        }
    }
for (size_t i = 0; i < stats.size(); i++) {
        center_stat[i] += stats[i]->centered_stat(point);
    }
}

std::unique_ptr<mvn_test> mvn_test_fixed::clone() {
    return std::make_unique<mvn_test_fixed>(*this);
}

mvn_test_fixed::mvn_test_fixed(const mvn_test_fixed& other) :mvn_test(other), stats(other.stats) {}

void mvn_test_fixed::compute_statistics() {}

mvn_test_gen::mvn_test_gen(std::shared_ptr<const Matrix> X, const Clustering& clst) : mvn_test(X, clst), X(X) {}

mvn_test_gen::mvn_test_gen(const mvn_test_gen& other) :mvn_test(other), X(other.X) {}

void mvn_test_gen::compute_statistics() {
    std::vector<unsigned> columns;
    for (size_t s: subset) {
        columns.insert(columns.end(), clustering->elements(s).begin(), clustering->elements(s).end());
    }
    std::shared_ptr<Matrix> Xs = std::make_shared<Matrix>(X->rows(), columns.size());
    for (int i = 0; i < columns.size(); i++) {
        Xs->col(i) = X->col(columns[i]);
    }
    Vector mean = Xs->rowwise().mean();
    Matrix centered = Xs->colwise() - mean;
    Matrix cov = (centered * centered.adjoint()) / double(centered.cols() - 1);
    mahalanobis_distances distances(Xs, cov, mean);

    for (unsigned ibeta = 0; ibeta < betas.size(); ibeta++) {
        double beta = betas[ibeta];
        mvn_stats stats(distances, *clustering, beta);
        center_stat[ibeta] = 0.0;
        auto n = Xs->cols();
        for (unsigned i = 0; i < n; i++) {
            center_stat[ibeta] += stats.centered_stat(i);
            for (unsigned j = 0; j < i; j++) {
                pairwise_stat[ibeta] += 2 *  stats.pairwise_stat(i, j);
            }
            pairwise_stat[ibeta] += stats.pairwise_stat(i, i);
        }
    }
}

void mvn_test_gen::remove(unsigned int i) {}

void mvn_test_gen::add(unsigned int i) {}

std::unique_ptr<mvn_test> mvn_test_gen::clone() {
    mvn_test_gen copy(*this);
    return std::make_unique<mvn_test_gen>(copy);
}
}