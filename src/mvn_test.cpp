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

    betas = {0.4, 0.8, 1.6, 3.2};

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

double mvn_test::get_normality_statistic() {
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
     mahalanobis_pairwise(Matrix::Zero(clst.size(), clst.size())) {
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

    std::vector<int> ids;
    for (int cluster: the_rest) {
        for (int el: clustering->elements(cluster)) {
            ids.push_back(el);
        }
    }

    auto lls = loglikelihood(ids);

    // find the cluster closest to the center
    std::vector<double> cl_lls;
    int i = 0;
    for (int cl: the_rest) {
        double mean = 0.0;
        for (int j = 0; j < clustering->cluster_size(cl); j++) {
            mean += lls[i + j];
        }
        mean /= (double)clustering->cluster_size(cl);
        i += clustering->cluster_size(cl);
        cl_lls.push_back(mean);
    }
    int pos = std::distance(cl_lls.begin(), std::max_element(cl_lls.begin(), cl_lls.end()));

    size_t point = the_rest[pos];
    effect_size += clustering->cluster_size(point);
    std::swap(the_rest[pos], the_rest[the_rest.size() - 1]);
    the_rest.pop_back();

    subset.push_back(point);
    add(point);
}

bool operator<(mvn_test& lhs, mvn_test& rhs) {
    return lhs.get_normality_statistic() < rhs.get_normality_statistic();
}

mvn_test::mvn_test(const mvn_test& other)
    :center_stat(other.center_stat),
     pairwise_stat(other.pairwise_stat),
     betas(other.betas),
     clustering(other.clustering),
     p(other.p),
     n(other.n),
     effect_size(other.effect_size),
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

const std::vector<int>& Clustering::elements(size_t i) const {
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
    distances = X->transpose() * S_inv * *X;
    mumu = mean.transpose() * S_inv * mean;
    diag = distances.diagonal();
}

double mahalanobis_distances::interpoint_distance(unsigned i, unsigned j) const {
    return diag(i) - 2 *  distances(i, j) + diag(j);
}

double mahalanobis_distances::distance(unsigned el) const {
    return diag(el) - ximu(el) - muxi(el) + mumu;
}

mvn_test_fixed::mvn_test_fixed(std::shared_ptr<const Matrix> X, const Clustering& clst, const Matrix& S, const Vector& mean)
        : mvn_test(X, clst), distances{std::make_shared<mahalanobis_distances>(X, S, mean)} {

    for (double beta: betas) {
        stats.push_back(std::make_shared<mvn_stats>(*distances, clst, beta));
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

mvn_test_fixed::mvn_test_fixed(const mvn_test_fixed& other) :mvn_test(other), distances(other.distances), stats(other.stats) {}

void mvn_test_fixed::compute_statistics() {}

std::vector<double> mvn_test_fixed::loglikelihood(const std::vector<int>& ids) const {
    std::vector<double> ret;
    for (int id: ids) {
        ret.push_back(-0.5 * distances->distance(id));
    }
    return ret;
}

mvn_test_gen::mvn_test_gen(std::shared_ptr<const Matrix> X, const Clustering& clst) : mvn_test(X, clst), X(X) {}

mvn_test_gen::mvn_test_gen(const mvn_test_gen& other) :mvn_test(other), X(other.X) {
}

namespace {
    std::vector<int> cluster_fold(std::shared_ptr<Clustering> clustering, const std::vector<size_t>& clusters) {
        std::vector<int> ret;
        for (int cluster: clusters) {
            auto& cls = clustering->elements(cluster);
            ret.insert(ret.end(), cls.begin(), cls.end());
        }
        return ret;
    }

    std::shared_ptr<Matrix> matrix_subset(std::shared_ptr<const Matrix> X, std::vector<int> columns) {
        std::shared_ptr<Matrix> Xs = std::make_shared<Matrix>(X->rows(), columns.size());
        for (int i = 0; i < columns.size(); i++) {
            Xs->col(i) = X->col(columns[i]);
        }
        return Xs;
    }

    mahalanobis_distances subset_distances(std::shared_ptr<const Matrix> X, const std::vector<int>& subset,
                                           const std::vector<int>& columns, bool euclidean = false) {
        auto Xs = matrix_subset(X, columns);
        auto Xsubset = matrix_subset(X, subset);
        Vector mean = Xsubset->rowwise().mean();
        Matrix centered = Xsubset->colwise() - mean;
        if (euclidean) {
            Matrix cov = Eigen::MatrixXd::Identity(centered.rows(), centered.rows());
            return {Xs, cov, mean};
        } else {
            Matrix cov = (centered * centered.adjoint()) / double(centered.cols() - 1);
            return {Xs, cov, mean};
        }
    }
}

void mvn_test_gen::compute_statistics() {
    auto columns = cluster_fold(clustering, subset);

    auto distances = subset_distances(X, columns, columns);
    for (unsigned ibeta = 0; ibeta < betas.size(); ibeta++) {
        double beta = betas[ibeta];
        mvn_stats stats(distances, *clustering, beta);
        center_stat[ibeta] = 0.0;
        pairwise_stat[ibeta] = 0.0;
        auto n = columns.size();
        for (unsigned i = 0; i < n; i++) {
            center_stat[ibeta] += stats.centered_stat(i);
            for (unsigned j = 0; j < i; j++) {
                pairwise_stat[ibeta] += 2 *  stats.pairwise_stat(i, j);
            }
            pairwise_stat[ibeta] += stats.pairwise_stat(i, i);
        }
    }
}

void mvn_test_gen::remove(unsigned int i) {
    if (subsample_size() > p + 1) {
        compute_statistics();
    }
}

void mvn_test_gen::add(unsigned int i) {
    if (subsample_size() > p + 1) {
        compute_statistics();
    }
}

std::unique_ptr<mvn_test> mvn_test_gen::clone() {
    mvn_test_gen copy(*this);
    return std::make_unique<mvn_test_gen>(copy);
}

std::vector<double> mvn_test_gen::loglikelihood(const std::vector<int>& ids) const {
    if (subset.size() == 0) {
        std::vector<double> ret;
        std::uniform_real_distribution<double> random;
        for (int i = 0; i < ids.size(); i++) {
            ret.push_back(random(wheel));
        }
        return ret;
    }
    auto sub = cluster_fold(clustering, subset);
    auto distances = subset_distances(X, sub, ids, true);
    std::vector<double> ret;
    for (int i = 0; i < ids.size(); i++) {
        ret.push_back(-0.5 * distances.distance(i));
    }
    return ret;
}
}