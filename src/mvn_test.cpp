#include <iostream>
#include <cmath>
#include <unordered_set>
#include <limits>
#include "include/third-party/cxxpool.h"

#include "include/mvn_test.h"

#undef NDEBUG
#include <assert.h>

namespace mvn {

namespace {

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

Matrix invert_matrix_or_throw(const Matrix& m) {
    Eigen::FullPivHouseholderQR<Matrix> qr(m);
    if (!qr.isInvertible()) {
        throw std::logic_error("Non-invertible matrix");
    }
    return qr.inverse();
}

// Returns A such that A.transpose() * A == spd (up to numerical tolerance).
Matrix sqrt_factor_of_spd(const Matrix& spd) {
    Eigen::LLT<Matrix> llt(spd);
    if (llt.info() == Eigen::Success) {
        // U^T U = spd
        return llt.matrixU();
    }

    // Fallback: eigen-decomposition (robust when LLT fails due to numerical issues)
    Matrix sym = 0.5 * (spd + spd.transpose());
    Eigen::SelfAdjointEigenSolver<Matrix> es(sym);
    if (es.info() != Eigen::Success) {
        throw std::logic_error("Unable to factorize covariance inverse");
    }
    Vector eval = es.eigenvalues();
    Matrix evec = es.eigenvectors();

    // Clamp tiny negatives from numerical noise
    for (int i = 0; i < eval.size(); i++) {
        if (eval(i) < 0 && eval(i) > -1e-12) {
            eval(i) = 0;
        }
        if (eval(i) < 0) {
            throw std::logic_error("Covariance inverse is not positive semidefinite");
        }
    }
    Vector sqrt_eval = eval.array().sqrt();
    return sqrt_eval.asDiagonal() * evec.transpose();
}

size_t choose_rff_dim(size_t p, size_t n_clusters) {
    // Heuristic: enough features to approximate kernel while keeping memory bounded.
    // The output influences memory as O(rff_dim * n_clusters).
    size_t dim = std::max<size_t>(128, 16 * p);
    dim = std::min<size_t>(512, dim);
    // When clusters are extremely many, keep dim smaller to reduce RAM.
    if (n_clusters > 200000) {
        dim = std::min<size_t>(dim, 256);
    }
    if (n_clusters > 1000000) {
        dim = std::min<size_t>(dim, 128);
    }
    return dim;
}

}

static void validate_cluster_id(size_t cluster_id, size_t n_clusters) {
    if (cluster_id >= n_clusters) {
        throw std::out_of_range("Cluster id is out of range");
    }
}

mvn_test_rff::mvn_test_rff(std::shared_ptr<const Matrix> X,
                           const Clustering& clst,
                           const Matrix& S,
                           const Vector& mean,
                           size_t requested_rff_dim)
        :clustering(std::make_shared<Clustering>(clst)),
         wheel(std::random_device()()) {
    if (X->cols() == 0 || X->rows() == 0) {
        throw std::invalid_argument("Matrix is empty");
    }

    betas = {0.8};

    n = clst.size();
    p = X->rows();
    if (n == 0) {
        throw std::invalid_argument("Clustering is empty");
    }
    if (p == 0) {
        throw std::invalid_argument("Matrix has zero rows");
    }
    if (mean.size() != static_cast<Eigen::Index>(p)) {
        throw std::invalid_argument("Mean vector dimension mismatch");
    }
    if (S.rows() != static_cast<Eigen::Index>(p) || S.cols() != static_cast<Eigen::Index>(p)) {
        throw std::invalid_argument("Covariance dimension mismatch");
    }

    pairwise_stat.resize(betas.size(), 0.0);
    center_stat.resize(betas.size(), 0.0);
    effect_size = 0;
    latest_subset_point = -1;

    if (n <= p) {
        throw std::logic_error("Too few points.");
    }

    rff_dim = requested_rff_dim == 0 ? choose_rff_dim(p, n) : requested_rff_dim;
    if (rff_dim == 0) {
        throw std::invalid_argument("RFF dimension must be positive");
    }

    // Precompute whitening and shared RFF parameters. Per-cluster values are computed lazily.
    auto params = std::make_shared<RFFParams>();
    params->X = std::move(X);
    params->A = sqrt_factor_of_spd(invert_matrix_or_throw(S));
    params->mean_whitened = params->A * mean;
    params->scale = std::sqrt(2.0 / static_cast<double>(rff_dim));

    // Shared RFF base frequencies and phases; beta is applied as a scalar multiplier.
    std::normal_distribution<double> normal(0.0, 1.0);
    std::uniform_real_distribution<double> unif(0.0, 2.0 * M_PI);
    params->W.resize((Eigen::Index)rff_dim, (Eigen::Index)p);
    params->b.resize((Eigen::Index)rff_dim);
    for (Eigen::Index i = 0; i < params->W.rows(); i++) {
        for (Eigen::Index j = 0; j < params->W.cols(); j++) {
            params->W(i, j) = normal(wheel);
        }
        params->b(i) = unif(wheel);
    }
    rff = std::move(params);

    subset_rff_sum.clear();
    subset_rff_sum.reserve(betas.size());
    for (size_t bi = 0; bi < betas.size(); bi++) {
        subset_rff_sum.emplace_back(Eigen::VectorXd::Zero((Eigen::Index)rff_dim));
    }

    cluster_rff_cache.clear();
    cluster_center_cache.clear();
    cluster_rff_cache.resize(betas.size());
    cluster_center_cache.resize(betas.size());

    // Uniform sampler over clusters (log-weights all zero).
    std::vector<double> lls(n, 0.0);
    sampler = RandomSampler(lls, wheel());

    while (effect_size < p + 1) {
        add_one();
    }
}

void mvn_test_rff::ensure_cluster_cached(size_t cluster_id) {
    validate_cluster_id(cluster_id, n);
    if (betas.empty()) {
        throw std::logic_error("No betas configured");
    }
    // If already cached for the first beta, assume cached for all.
    if (cluster_rff_cache[0].find(cluster_id) != cluster_rff_cache[0].end()) {
        return;
    }

    std::vector<Eigen::VectorXf> phi_sums;
    std::vector<float> center_sums;
    phi_sums.reserve(betas.size());
    center_sums.reserve(betas.size());
    for (size_t bi = 0; bi < betas.size(); bi++) {
        phi_sums.emplace_back(Eigen::VectorXf::Zero((Eigen::Index)rff_dim));
        center_sums.emplace_back(0.0f);
    }

    for (int el : clustering->elements(cluster_id)) {
        if (el < 0 || el >= rff->X->cols()) {
            throw std::out_of_range("Cluster element index is out of range");
        }
        Vector y = rff->A * rff->X->col(el);
        Vector proj = rff->W * y; // rff_dim
        Vector y_centered = y - rff->mean_whitened;
        const double d2_to_mean = y_centered.squaredNorm();

        for (size_t bi = 0; bi < betas.size(); bi++) {
            const double beta = betas[bi];
            const double k_center = -(beta * beta / (2.0 * (1.0 + beta * beta)));
            center_sums[bi] += (float)std::exp(k_center * d2_to_mean);

            auto& phi = phi_sums[bi];
            for (Eigen::Index d = 0; d < (Eigen::Index)rff_dim; d++) {
                phi(d) += (float)(rff->scale * std::cos(beta * proj(d) + rff->b(d)));
            }
        }
    }

    for (size_t bi = 0; bi < betas.size(); bi++) {
        cluster_rff_cache[bi].emplace(cluster_id, std::move(phi_sums[bi]));
        cluster_center_cache[bi].emplace(cluster_id, center_sums[bi]);
    }
}

double mvn_test_rff::get_normality_statistic() {
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

const std::vector<size_t>& mvn_test_rff::current_subset() const {
    return subset;
}

size_t mvn_test_rff::dimensions() const {
    return p;
}

size_t mvn_test_rff::subsample_size() const {
    return subset.size();
}

size_t mvn_test_rff::sample_size() const {
    return n;
}

void mvn_test_rff::swap_once(bool reject_last) {
    if (subset.empty() || sampler.n_active() == 0) {
        throw std::logic_error("Unable to swap points.");
    }

    int replacing_point = -1;
    if (reject_last) {
        assert(latest_subset_point != -1);
        replacing_point = latest_subset_point;
        latest_subset_point = -1;
    } else {
        std::uniform_int_distribution<unsigned> subset_unif(0, subsample_size() - 1);
        std::swap(subset[subset_unif(wheel)], subset.back());
        replacing_point = sampler.sample();
    }

    auto subset_point = subset.back();
    if (!reject_last) {
        latest_subset_point = subset_point;
    }

    effect_size += clustering->cluster_size(replacing_point) - clustering->cluster_size(subset_point);

    remove(subset_point);

    sampler.enable(subset.back());
    sampler.disable(replacing_point);

    subset.pop_back();
    subset.push_back(replacing_point);

    add(replacing_point);
}

void mvn_test_rff::add_one() {
    if (sampler.n_active() == 0) {
        throw std::logic_error("Can't add point to the model.");
    }
    latest_subset_point = -1;

    size_t point = sampler.sample();
    effect_size += clustering->cluster_size(point);
    sampler.disable(point);

    subset.push_back(point);
    add(point);
}

mvn_test_rff::mvn_test_rff(const mvn_test_rff& other)
    :sampler(other.sampler),
     pairwise_stat(other.pairwise_stat),
     center_stat(other.center_stat),
     betas(other.betas),
     rff_dim(other.rff_dim),
     rff(other.rff),
     subset_rff_sum(other.subset_rff_sum),
     cluster_rff_cache(other.cluster_rff_cache),
     cluster_center_cache(other.cluster_center_cache),
     clustering(other.clustering),
     p(other.p),
     n(other.n),
     effect_size(other.effect_size),
     latest_subset_point(other.latest_subset_point),
     wheel{other.wheel()},
     subset(other.subset) {}

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

void mvn_test_rff::remove(unsigned point) {
    ensure_cluster_cached(point);
    for (size_t bi = 0; bi < betas.size(); bi++) {
        subset_rff_sum[bi] -= cluster_rff_cache[bi].at(point).cast<double>();
        pairwise_stat[bi] = subset_rff_sum[bi].squaredNorm();
        center_stat[bi] -= (double)cluster_center_cache[bi].at(point);
    }
}

void mvn_test_rff::add(unsigned int point) {
    ensure_cluster_cached(point);
    for (size_t bi = 0; bi < betas.size(); bi++) {
        subset_rff_sum[bi] += cluster_rff_cache[bi].at(point).cast<double>();
        pairwise_stat[bi] = subset_rff_sum[bi].squaredNorm();
        center_stat[bi] += (double)cluster_center_cache[bi].at(point);
    }
}

std::shared_ptr<mvn_test_base> mvn_test_rff::clone() const {
    return std::make_shared<mvn_test_rff>(*this);
}

RandomSampler::RandomSampler(const std::vector<double>& logscale, long seed) :runif(0.0, 1.0), wheel(seed), original(logscale),
                                                                              size(logscale.size()) {
    if (original.size() < 2) {
        throw std::invalid_argument("Too little segment tree");
    }

    active_tree.resize(2 * original.size() - 1);
    segment_tree.resize(2 * original.size() - 1);
    for (size_t i = 0; i < original.size(); i++) {
        auto k = segment_tree.size() - i - 1;
        active_tree[k] = 1;
        segment_tree[k] = original[i];
    }
    for (int i = segment_tree.size() - size - 1; i >= 0; i--) {
        update_inner_node(i);
    }
}


bool RandomSampler::is_active(size_t n) const {
    return active_tree.at(el_pos(n));
}

void RandomSampler::disable(size_t n) {
    if (!is_active(n)) {
        throw std::logic_error("Disabling non-active element");
    }
    auto pos = el_pos(n);
    segment_tree[pos] = -std::numeric_limits<double>::infinity();
    active_tree[pos] = false;
    update(pos);
}

void RandomSampler::enable(size_t n) {
    if (is_active(n)) {
        throw std::logic_error("Enabling active element.");
    }
    auto pos = el_pos(n);
    segment_tree[pos] = original[n];
    active_tree.at(pos) = true;
    update(pos);
}

size_t RandomSampler::sample() {
    if (size < 2) {
        throw std::logic_error("Too little segment tree");
    }
    size_t node = 0;
    while (!is_leaf(node)) {
        auto chld = children(node);
        if (active_tree.at(chld.first) == 0) {
            if (active_tree.at(chld.second) == 0) {
                throw std::logic_error("No active elements in subtree");
            }
            node = chld.second;
        } else if(active_tree.at(chld.second) == 0) {
            node = chld.first;
        } else {
            double l = segment_tree.at(chld.first);
            double r = segment_tree.at(chld.second);
            double maxL = std::max(l, r);
            l = std::exp(l - maxL);
            r = std::exp(r - maxL);
            auto sum = l + r;
            l /= sum;
            if (runif(wheel) < l) {
                node = chld.first;
            } else {
                node = chld.second;
            }
        }
    }

    int aug_nodes = segment_tree.size() - original.size();
    if (!active_tree.at(node)) {
        throw std::logic_error("Sampled element is not active.");
    }
    int element = node - aug_nodes;
    assert(el_pos(element) == node);
    return element;
}

std::pair<size_t, size_t> RandomSampler::children(size_t node) const {
    if (node > segment_tree.size() - original.size() - 1) {
        throw std::invalid_argument("It's a leaf");
    }
    return {2 * node + 1, 2 * node + 2};
}

bool RandomSampler::is_root(size_t node) {
    return node == 0;
}

bool RandomSampler::is_leaf(size_t node) const {
    return node >= segment_tree.size() - original.size();
}

size_t RandomSampler::el_pos(size_t el) const {
    size_t pos = segment_tree.size() -  original.size() + el;
    if (pos >= segment_tree.size() || pos < 0) {
        throw std::out_of_range("Out of range");
    }
    return pos;
}

double RandomSampler::sum_log(double l, double r) {
    if (std::isinf(l)) {
        return r;
    }
    if (std::isinf(r)) {
        return l;
    }
    double maxL = std::max(l, r);
    l = std::exp(l - maxL);
    r = std::exp(r - maxL);
    return std::log(l + r) + maxL;
}

size_t RandomSampler::parent(size_t node) {
    return (node - 1) / 2;
}

void RandomSampler::update_inner_node(size_t node) {
    auto chs = children(node);
    segment_tree[node] = sum_log(segment_tree[chs.first], segment_tree[chs.second]);
    assert(!std::isnan(segment_tree[node]));
    active_tree[node] = active_tree.at(chs.first) + active_tree.at(chs.second);
}

void RandomSampler::update(size_t node) {
    if (node < segment_tree.size() - original.size()) {
        throw std::invalid_argument("That's not a leaf");
    }
    node = parent(node);
    while(true) {
        update_inner_node(node);
        if (is_root(node)) {
            break;
        }
        node = parent(node);
    }
}

RandomSampler::RandomSampler() :runif(0.0, 1.0), wheel(0), size(0) {}

RandomSampler::RandomSampler(const RandomSampler& other) :runif(0.0, 1.0), wheel(other.wheel()),
                                                          original(other.original), segment_tree(other.segment_tree),
                                                          active_tree(other.active_tree), size(other.size) {}

RandomSampler& RandomSampler::operator=(RandomSampler&& other) {
    runif = other.runif;
    wheel = other.wheel;
    original = std::move(other.original);
    segment_tree = std::move(other.segment_tree);
    active_tree = std::move(other.active_tree);
    size = other.size;
    return *this;
}

size_t RandomSampler::n_active() const {
    return active_tree.at(0);
}

}
