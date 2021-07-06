#include <iostream>
#include <cmath>
#include <unordered_set>

#include "include/mvn_test.h"

namespace mvn {

mvn_test::mvn_test(std::shared_ptr<const Matrix> X, const Clustering& clst, const Matrix& S, const Vector& mean)
        :distances{std::make_shared<mahalanobis_distances>(X, S, mean)},
         clustering(std::make_shared<Clustering>(clst)),
         wheel(std::random_device()()) {
    if (X->cols() == 0 || X->rows() == 0) {
        throw std::invalid_argument("Matrix is empty");
    }

    betas = {0.4, 0.8, 1.6, 3.2};

    pairwise_stat.resize(betas.size(), 0.0);
    center_stat.resize(betas.size(), 0.0);

    n = clst.size();
    p = X->rows();
    effect_size = 0;
    latest_subset_point = -1;

    if (n <= p) {
        throw std::logic_error("Too few points.");
    }

    for (double beta: betas) {
        stats.push_back(std::make_shared<mvn_stats>(*distances, clst, beta));
    }

    while (effect_size < p + 1) {
        add_one();
    }

    std::vector<double> lls;
    for (size_t i = 0; i < n; i++) {
        auto ids = clst.elements(i);
        auto loglikelihoods = loglikelihood(ids);
        lls.push_back(std::accumulate(loglikelihoods.begin(), loglikelihoods.end(), 0.0));
    }

    sampler = RandomSampler(lls, wheel());
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
        latest_subset_point = replacing_point;
    }

    auto subset_point = subset.back();

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

void mvn_test::remove(unsigned point) {
    for (size_t i = 0; i < stats.size(); i++) {
        for (auto s: subset) {
            pairwise_stat[i] -= 2 * stats[i]->pairwise_stat(point, s);
        }
    }

    for (size_t i = 0; i < stats.size(); i++) {
        center_stat[i] -= stats[i]->centered_stat(point);
    }
}

void mvn_test::add(unsigned int point) {
    for (size_t i = 0; i < stats.size(); i++) {
        for (auto s: subset) {
            pairwise_stat[i] += 2 * stats[i]->pairwise_stat(point, s);
        }
    }
    for (size_t i = 0; i < stats.size(); i++) {
        center_stat[i] += stats[i]->centered_stat(point);
    }
}

std::unique_ptr<mvn_test> mvn_test::clone() {
    return std::make_unique<mvn_test>(*this);
}

std::vector<double> mvn_test::loglikelihood(const std::vector<int>& ids) const {
    std::vector<double> ret;
    for (int id: ids) {
        ret.push_back(-0.5 * distances->distance(id));
    }
    return ret;
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
    return node - aug_nodes;
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
    size_t pos = segment_tree.size() -  original.size() - 1 + el;
    if (pos >= segment_tree.size()) {
        throw std::out_of_range("Out of range");
    }
    return pos;
}

double RandomSampler::sum_log(double l, double r) {
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
    runif = std::move(other.runif);
    wheel = std::move(other.wheel);
    original = std::move(other.original);
    segment_tree = std::move(other.segment_tree);
    active_tree = std::move(other.active_tree);
    return *this;
}

size_t RandomSampler::n_active() const {
    return active_tree.at(0);
}

}
