#include <algorithm>
#include <iostream>
#include <cmath>
#include <unordered_set>
#include "include/third-party/cxxpool.h"

#include "include/mvn_test.h"

#undef NDEBUG
#include <assert.h>

namespace mvn {

namespace {

size_t precompute_workers(size_t tasks, size_t requested_threads) {
    if (tasks == 0) {
        return 1;
    }

    size_t workers = requested_threads;
    if (workers == 0) {
        workers = std::thread::hardware_concurrency();
    }
    if (workers == 0) {
        workers = 1;
    }

    return std::min(tasks, workers);
}

size_t effective_cluster_tile_size(size_t requested_tile_size) {
    return requested_tile_size == 0 ? 32 : requested_tile_size;
}

size_t packed_pairwise_index(size_t i, size_t j, size_t n) {
    if (i > j) {
        std::swap(i, j);
    }
    return i * n - (i * (i - 1)) / 2 + (j - i);
}

Matrix gather_columns(const Matrix& matrix, const std::vector<int>& cluster, size_t begin, size_t end) {
    Matrix gathered(matrix.rows(), end - begin);
    for (size_t offset = 0; offset < end - begin; offset++) {
        gathered.col(offset) = matrix.col(cluster[begin + offset]);
    }

    return gathered;
}

Eigen::ArrayXd gather_quadratic_forms(const std::vector<double>& quadratic_forms,
                                      const std::vector<int>& cluster,
                                      size_t begin,
                                      size_t end) {
    Eigen::ArrayXd gathered(end - begin);
    for (size_t offset = 0; offset < end - begin; offset++) {
        gathered(offset) = quadratic_forms[cluster[begin + offset]];
    }

    return gathered;
}

double accumulate_centered_cluster(const mahalanobis_distances& distances,
                                   const std::vector<int>& cluster,
                                   double k_center) {
    double centered = 0.0;
    for (int el: cluster) {
        centered += std::exp(k_center * distances.distance(el));
    }

    return centered;
}

double accumulate_pairwise_cluster(const mahalanobis_distances& distances,
                                   const std::vector<int>& left,
                                   const std::vector<int>& right,
                                   double k_pw,
                                   bool same_cluster,
                                   size_t cluster_tile_size) {
    double sum = 0.0;

    const auto& samples = distances.samples();
    const auto& transformed = distances.transformed_samples();
    const auto& quadratic_forms = distances.quadratic_forms();

    if (same_cluster) {
        for (size_t left_begin = 0; left_begin < left.size(); left_begin += cluster_tile_size) {
            size_t left_end = std::min(left.size(), left_begin + cluster_tile_size);
            Matrix left_samples = gather_columns(samples, left, left_begin, left_end);
            Eigen::ArrayXd left_quad = gather_quadratic_forms(quadratic_forms, left, left_begin, left_end);

            for (size_t right_begin = left_begin; right_begin < left.size(); right_begin += cluster_tile_size) {
                size_t right_end = std::min(left.size(), right_begin + cluster_tile_size);
                Matrix right_transformed = gather_columns(transformed, left, right_begin, right_end);
                Eigen::ArrayXd right_quad = gather_quadratic_forms(quadratic_forms, left, right_begin, right_end);

                Matrix cross(left_end - left_begin, right_end - right_begin);
                cross.noalias() = left_samples.transpose() * right_transformed;
                cross *= -2.0;
                cross.colwise() += left_quad.matrix();
                cross.rowwise() += right_quad.transpose().matrix();
                Eigen::ArrayXXd kernel = cross.array();
                kernel *= k_pw;
                kernel = kernel.exp();

                if (left_begin == right_begin) {
                    for (Eigen::Index i = 0; i < kernel.rows(); i++) {
                        sum += 0.5 * kernel(i, i);
                        for (Eigen::Index j = i + 1; j < kernel.cols(); j++) {
                            sum += kernel(i, j);
                        }
                    }
                } else {
                    sum += kernel.sum();
                }
            }
        }

        return sum;
    }

    for (size_t left_begin = 0; left_begin < left.size(); left_begin += cluster_tile_size) {
        size_t left_end = std::min(left.size(), left_begin + cluster_tile_size);
        Matrix left_samples = gather_columns(samples, left, left_begin, left_end);
        Eigen::ArrayXd left_quad = gather_quadratic_forms(quadratic_forms, left, left_begin, left_end);

        for (size_t right_begin = 0; right_begin < right.size(); right_begin += cluster_tile_size) {
            size_t right_end = std::min(right.size(), right_begin + cluster_tile_size);
            Matrix right_transformed = gather_columns(transformed, right, right_begin, right_end);
            Eigen::ArrayXd right_quad = gather_quadratic_forms(quadratic_forms, right, right_begin, right_end);

            Matrix cross(left_end - left_begin, right_end - right_begin);
            cross.noalias() = left_samples.transpose() * right_transformed;
            cross *= -2.0;
            cross.colwise() += left_quad.matrix();
            cross.rowwise() += right_quad.transpose().matrix();
            Eigen::ArrayXXd kernel = cross.array();
            kernel *= k_pw;
            sum += kernel.exp().sum();
        }
    }

    return sum;
}

void precompute_cluster_range(std::vector<double>& mahalanobis_centered,
                              std::vector<double>& mahalanobis_pairwise,
                              size_t n_clusters,
                              const mahalanobis_distances& distances,
                              const Clustering& clst,
                              double k_center,
                              double k_pw,
                              size_t begin,
                              size_t end,
                              size_t cluster_tile_size) {
    for (size_t cl = begin; cl < end; cl++) {
        const auto& cluster = clst.elements(cl);
        mahalanobis_centered[cl] = accumulate_centered_cluster(distances, cluster, k_center);
        for (size_t pair_cl = cl; pair_cl < n_clusters; pair_cl++) {
            double value = accumulate_pairwise_cluster(distances,
                                                       cluster,
                                                       clst.elements(pair_cl),
                                                       k_pw,
                                                       cl == pair_cl,
                                                       cluster_tile_size);
            mahalanobis_pairwise[packed_pairwise_index(cl, pair_cl, n_clusters)] = value;
        }
    }
}

}

mvn_test::mvn_test(std::shared_ptr<const Matrix> X, const Clustering& clst, const Matrix& S, const Vector& mean,
           const PrecomputeConfig& config, std::mt19937::result_type seed)
        :distances{std::make_shared<mahalanobis_distances>(X, S, mean)},
         clustering(std::make_shared<Clustering>(clst)),
         wheel(seed) {
    if (X->cols() == 0 || X->rows() == 0) {
        throw std::invalid_argument("Matrix is empty");
    }

    betas = {0.8};

    pairwise_stat.resize(betas.size(), 0.0);
    center_stat.resize(betas.size(), 0.0);

    n = clst.size();
    p = X->rows();
    effect_size = 0;
    latest_subset_point = -1;

    if (n <= p) {
        throw std::logic_error("Too few points.");
    }

    stats.reserve(betas.size());
    for (double beta: betas) {
        stats.push_back(std::make_shared<mvn_stats>(*distances, clst, beta, config));
    }

    std::vector<double> lls(n, 0.0);

    sampler = RandomSampler(lls, wheel());
    distances.reset();

    while (effect_size < p + 1) {
        add_one();
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

mvn_stats::mvn_stats(const mahalanobis_distances& distances, const Clustering& clst, double beta,
                     const PrecomputeConfig& config)
    :mahalanobis_centered(clst.size(), 0.0),
     mahalanobis_pairwise(clst.size() * (clst.size() + 1) / 2, 0.0),
     n_clusters(clst.size()) {
    size_t n = clst.size();
    size_t cluster_tile_size = effective_cluster_tile_size(config.cluster_tile_size);

    double k_center = -(beta * beta / (2 * (1 + beta * beta)));
    double k_pw = -(beta * beta) / 2.0;

    size_t workers = precompute_workers(n, config.threads);
    if (workers == 1) {
        precompute_cluster_range(mahalanobis_centered,
                                 mahalanobis_pairwise,
                                 n,
                                 distances,
                                 clst,
                                 k_center,
                                 k_pw,
                                 0,
                                 n,
                                 cluster_tile_size);
        return;
    }

    size_t block_size = (n + workers - 1) / workers;

    cxxpool::thread_pool pool(workers);
    std::vector<std::future<void>> futures;
    futures.reserve(workers);
    for (size_t block_begin = 0; block_begin < n; block_begin += block_size) {
        size_t block_end = std::min(n, block_begin + block_size);
        futures.push_back(pool.push([this, &clst, &distances, k_center, k_pw, block_begin, block_end, cluster_tile_size, n]() {
            precompute_cluster_range(mahalanobis_centered,
                                     mahalanobis_pairwise,
                                     n,
                                     distances,
                                     clst,
                                     k_center,
                                     k_pw,
                                     block_begin,
                                     block_end,
                                     cluster_tile_size);
        }));
    }

    cxxpool::get(futures.begin(), futures.end());
}

double mvn_stats::pairwise_stat(size_t i, size_t j) const {
    return mahalanobis_pairwise[packed_pairwise_index(i, j, n_clusters)];
}

double mvn_stats::centered_stat(size_t i) const {
    return mahalanobis_centered[i];
}

double mvn_stats::sum_pairwise(size_t point, const std::vector<size_t>& ss) const {
    double ret = 0.0;
    for (auto s: ss) {
        ret += mahalanobis_pairwise[packed_pairwise_index(point, s, n_clusters)];
    }
    return ret;
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

void mvn_test::add_one() {
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

bool operator<(mvn_test& lhs, mvn_test& rhs) {
    return lhs.get_normality_statistic() < rhs.get_normality_statistic();
}

mvn_test::mvn_test(const mvn_test& other)
    :distances(other.distances),
     stats(other.stats),
     sampler(other.sampler),
     pairwise_stat(other.pairwise_stat),
     center_stat(other.center_stat),
     betas(other.betas),
     clustering(other.clustering),
     p(other.p),
     n(other.n),
     effect_size(other.effect_size),
     latest_subset_point(other.latest_subset_point),
     wheel{other.wheel},
     subset(other.subset) {}

void mvn_test::reseed(std::mt19937::result_type seed) {
    wheel.seed(seed);
    sampler.reseed(wheel());
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

const std::vector<int>& Clustering::elements(size_t i) const {
    return clusters.at(i);
}

size_t Clustering::cluster_size(size_t i) const {
    return clusters.at(i).size();
}

mahalanobis_distances::mahalanobis_distances(std::shared_ptr<const Matrix> X, const Matrix& S, const Vector& mean)
        :X(std::move(X)),
         transformed(S.rows(), this->X->cols()),
         quadratic_form(this->X->cols()),
         centered_distance(this->X->cols()) {
    if (S.rows() != S.cols() || S.rows() != this->X->rows() || mean.size() != this->X->rows()) {
        throw std::invalid_argument("Covariance and mean dimensions must match the number of components (rows) of X.");
    }
    Eigen::LDLT<Matrix> solver(S);
    if (solver.info() != Eigen::Success) {
        throw std::logic_error("Non-invertible matrix. Must not happen.");
    }

    transformed = solver.solve(*this->X);
    if (solver.info() != Eigen::Success) {
        throw std::logic_error("Failed to solve covariance system. Must not happen.");
    }

    Vector transformed_mean = solver.solve(mean);
    if (solver.info() != Eigen::Success) {
        throw std::logic_error("Failed to solve centered covariance system. Must not happen.");
    }
    double mumu = mean.dot(transformed_mean);

    for (Eigen::Index i = 0; i < this->X->cols(); i++) {
        quadratic_form[i] = this->X->col(i).dot(transformed.col(i));
        double ximu = mean.dot(transformed.col(i));
        centered_distance[i] = quadratic_form[i] - 2 * ximu + mumu;
    }
}

double mahalanobis_distances::interpoint_distance(unsigned i, unsigned j) const {
    return quadratic_form[i] - 2 * X->col(i).dot(transformed.col(j)) + quadratic_form[j];
}

double mahalanobis_distances::distance(unsigned el) const {
    return centered_distance[el];
}

const Matrix& mahalanobis_distances::samples() const {
    return *X;
}

const Matrix& mahalanobis_distances::transformed_samples() const {
    return transformed;
}

const std::vector<double>& mahalanobis_distances::quadratic_forms() const {
    return quadratic_form;
}

void mvn_test::remove(unsigned point) {
    for (size_t i = 0; i < stats.size(); i++) {
        pairwise_stat[i] -= 2 * stats[i]->sum_pairwise(point, subset);
    }

    for (size_t i = 0; i < stats.size(); i++) {
        center_stat[i] -= stats[i]->centered_stat(point);
    }
}

void mvn_test::add(unsigned int point) {
    for (size_t i = 0; i < stats.size(); i++) {
        pairwise_stat[i] += 2 * stats[i]->sum_pairwise(point, subset);
    }
    for (size_t i = 0; i < stats.size(); i++) {
        center_stat[i] += stats[i]->centered_stat(point);
    }
}

std::unique_ptr<mvn_test> mvn_test::clone(std::mt19937::result_type seed) const {
    auto copy = std::make_unique<mvn_test>(*this);
    copy->reseed(seed);
    return copy;
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

RandomSampler::RandomSampler(const RandomSampler& other) :runif(other.runif), wheel(other.wheel),
                                                          original(other.original), segment_tree(other.segment_tree),
                                                          active_tree(other.active_tree), size(other.size) {}

void RandomSampler::reseed(std::mt19937::result_type seed) {
    wheel.seed(seed);
}

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
