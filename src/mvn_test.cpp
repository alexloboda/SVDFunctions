#include <iostream>
#include <cmath>
#include <unordered_set>
#include "include/third-party/cxxpool.h"
#include "include/exceptions.h"
#include "include/RandomSampler.h"
#include "include/mvn_clst.h"
#include "include/mahalanobis_distances.h"
#include "include/mvn_test.h"
#include "include/mvn_stats.h"
#include "include/mvn_stats_approx.h"
#include "include/mvn_stats_interpoint.h"

#undef NDEBUG
#include <assert.h>

namespace mvn {

mvn_test::mvn_test(std::shared_ptr<const Matrix> X, const Clustering& clst, const Matrix& S, const Vector& mean)
        :distances{std::make_shared<mahalanobis_distances>(X, S, mean)},
         clustering(std::make_shared<Clustering>(clst)),
         wheel(std::random_device()()) {
    initialize_common(X, clst, S, mean, std::nullopt);
}

mvn_test::mvn_test(std::shared_ptr<const Matrix> X, const Clustering& clst, const Matrix& S, const Vector& mean, const std::string& filename)
        :distances{std::make_shared<mahalanobis_distances>(X, S, mean)},
         clustering(std::make_shared<Clustering>(clst)),
         wheel(std::random_device()()) {
    initialize_common(X, clst, S, mean, filename);
}

void mvn_test::initialize_common(std::shared_ptr<const Matrix> X, const Clustering& clst, const Matrix& S, const Vector& mean, std::optional<std::string> filename) {
    if (X->cols() == 0 || X->rows() == 0) {
        throw std::invalid_argument("Matrix is empty");
    }

    betas = {0.3};
    pairwise_stat.resize(betas.size(), 0.0);
    center_stat.resize(betas.size(), 0.0);

    n = clst.size();
    p = X->rows();
    effect_size = 0;
    latest_subset_point = -1;

    if (n <= p) {
        throw std::logic_error("Too few points.");
    }

    initialize_stats(clst, filename);

    std::vector<double> lls = compute_loglikelihoods(clst);

    sampler = RandomSampler(lls, wheel());

    while (effect_size < p + 1) {
        add_one();
    }
}

void mvn_test::initialize_stats(const Clustering& clst, std::optional<std::string> filename) {
    for (double beta : betas) {
        std::shared_ptr<mvn_stats> stat = filename.has_value()
            ? std::static_pointer_cast<mvn_stats>(std::make_shared<mvn_stats_approx>(filename.value(), n))
            : std::static_pointer_cast<mvn_stats>(std::make_shared<mvn_stats_interpoint>(n));
        stat->init_pairwise(*distances, clst, beta);
        stats.push_back(stat);
    }
}

std::vector<double> mvn_test::compute_loglikelihoods(const Clustering& clst) const {
    std::vector<double> lls;
    for (size_t i = 0; i < n; i++) {
        auto ids = clst.elements(i);
        auto loglikelihoods = loglikelihood(ids);
        lls.push_back(0.0); // Placeholder for actual computation
    }
    return lls;
}

void mvn_test::update_stats(unsigned point, bool is_addition) {
    int factor = is_addition ? 1 : -1;
    for (size_t i = 0; i < stats.size(); i++) {
        pairwise_stat[i] += factor * 2 * stats[i]->sum_pairwise(point, subset);
        center_stat[i] += factor * stats[i]->centered_stat(point);
    }
}

void mvn_test::remove(unsigned point) {
    update_stats(point, false);
}

void mvn_test::add(unsigned int point) {
    update_stats(point, true);
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
     wheel{other.wheel()},
     subset(other.subset) {}

std::unique_ptr<mvn_test> mvn_test::clone() {
    return std::make_unique<mvn_test>(*this);
}

std::vector<double> mvn_test::loglikelihood(const std::vector<int>& ids) const {
    std::vector<double> ret;
    for (int id: ids) {
        double dist = std::sqrt(distances->distance(id));
        if (dist > 1) {
            dist = 1.0 / (dist * dist);
        } else {
            dist = 1;
        }
        ret.push_back(std::log(dist));
    }
    return ret;
}

}
