#include <iostream>
#include <cmath>
#include <unordered_set>
#include <numeric>
#include <algorithm>
#include "include/third-party/cxxpool.h"

#include "include/mvn_test.h"

#if defined(__SSE__)
#include <xmmintrin.h>
#endif

namespace mvn {

namespace {

constexpr size_t NYSTROM_DEFAULT_MAX_FEATURES = 1024;
constexpr double NYSTROM_REL_EIGEN_FLOOR = 1e-6;
constexpr double NYSTROM_MAX_COND = 1e8;
constexpr double NYSTROM_RIDGE_SCALE = 1e-8;
constexpr double NYSTROM_RIDGE_MIN = 1e-10;

double dot_float_scalar(const float* lhs, const float* rhs, size_t n) {
    double acc = 0.0;
    for (size_t i = 0; i < n; ++i) {
        acc += static_cast<double>(lhs[i]) * static_cast<double>(rhs[i]);
    }
    return acc;
}

double dot_float_sse(const float* lhs, const float* rhs, size_t n) {
#if defined(__SSE__)
    size_t i = 0;
    __m128 acc = _mm_setzero_ps();
    for (; i + 4 <= n; i += 4) {
        __m128 a = _mm_loadu_ps(lhs + i);
        __m128 b = _mm_loadu_ps(rhs + i);
        acc = _mm_add_ps(acc, _mm_mul_ps(a, b));
    }
    alignas(16) float tmp[4];
    _mm_store_ps(tmp, acc);
    double sum = static_cast<double>(tmp[0]) + static_cast<double>(tmp[1]) +
                 static_cast<double>(tmp[2]) + static_cast<double>(tmp[3]);
    for (; i < n; ++i) {
        sum += static_cast<double>(lhs[i]) * static_cast<double>(rhs[i]);
    }
    return sum;
#else
    return dot_float_scalar(lhs, rhs, n);
#endif
}

std::vector<size_t> choose_landmarks(size_t n_samples, size_t m, std::mt19937& wheel) {
    std::vector<size_t> all(n_samples);
    std::iota(all.begin(), all.end(), 0);
    std::shuffle(all.begin(), all.end(), wheel);
    all.resize(m);
    return all;
}

std::vector<size_t> choose_rff_levels(size_t max_features) {
    std::vector<size_t> levels;
    if (max_features == 0) {
        return levels;
    }

    size_t level = std::min<size_t>(64, max_features);
    while (level < max_features) {
        levels.push_back(level);
        const size_t next = std::min(max_features, level * 2);
        if (next == level) {
            break;
        }
        level = next;
    }
    levels.push_back(max_features);
    return levels;
}

} // namespace

void mvn_test::check_aux_state() const {
    if (!use_hybrid) {
        return;
    }
    if (!rff_stats) {
        throw std::logic_error("Hybrid mode is enabled but RFF stats are missing");
    }
    if (rff_feature_levels.empty()) {
        throw std::logic_error("Hybrid mode is enabled but the RFF ladder is empty");
    }
    if (rff_pairwise_stats.size() != rff_feature_levels.size()) {
        throw std::logic_error("RFF ladder state is inconsistent: pairwise stats and levels differ in size");
    }
    if (rff_subset_feature_sum.size() != rff_stats->features_dim()) {
        throw std::logic_error("RFF ladder state is inconsistent: subset feature sum has unexpected size");
    }
}

mvn_test::mvn_test(std::shared_ptr<const Matrix> X, const Clustering& clst, const Matrix& S, const Vector& mean,
           bool use_nystrom, size_t n_features,
           bool use_hybrid, size_t rff_features, uint32_t seed)
        :distances{std::make_shared<mahalanobis_distances>(X, S, mean)},
         clustering(std::make_shared<Clustering>(clst)),
     use_nystrom(use_nystrom),
     n_features(n_features),
     use_hybrid(use_hybrid),
         wheel(std::random_device()()) {
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
    latest_replacing_point = -1;

    if (n <= p) {
        throw std::logic_error("Too few points.");
    }

    cxxpool::thread_pool pool(std::thread::hardware_concurrency());
    std::vector<std::future<std::shared_ptr<mvn_stats>>> futures;
    for (double beta: betas) {
        futures.push_back(pool.push([this, clustering = this->clustering, beta, seed]() -> std::shared_ptr<mvn_stats> {
            const uint32_t local_seed = seed ^ (uint32_t)(beta * 1e6);
            return std::make_shared<mvn_stats>(*distances, *clustering, beta, this->use_nystrom, this->n_features,
                                               false, 0, local_seed);
        }));
    }
    for (int i = 0; i < betas.size(); i++) {
        stats.push_back(futures[i].get());
    }

    pairwise_stat.assign(betas.size(), 0.0);
    center_stat.assign(betas.size(), 0.0);
    subset_feature_sum.clear();
    subset_feature_sum.resize(betas.size());
    for (size_t i = 0; i < stats.size(); ++i) {
        if (stats[i]->is_feature_mode()) {
            subset_feature_sum[i].assign(stats[i]->features_dim(), 0.0f);
        }
    }

    if (use_hybrid) {
        // Single-beta auxiliary statistic for gating with progressively larger
        // RFF prefixes taken from one precomputed embedding bank.
        const double beta = betas.front();
        const uint32_t rff_seed = seed + 1337u;
        rff_stats = std::make_shared<mvn_stats>(*distances, *this->clustering, beta,
                                                false, 0,
                                                true, rff_features, rff_seed);
        rff_feature_levels = choose_rff_levels(rff_stats->features_dim());
        rff_pairwise_stats.assign(rff_feature_levels.size(), 0.0);
        rff_center_stat = 0.0;
        if (rff_stats->is_feature_mode()) {
            rff_subset_feature_sum.assign(rff_stats->features_dim(), 0.0f);
        }
        check_aux_state();
    }

    std::vector<double> lls;
    for (size_t i = 0; i < n; i++) {
        auto ids = clst.elements(i);
        auto loglikelihoods = loglikelihood(ids);
        lls.push_back(0.0);
    }

    sampler = RandomSampler(lls, wheel());

    while (effect_size < p + 1) {
        add_one();
    }
}

double mvn_test::get_aux_normality_statistic() {
    if (!has_aux_statistic()) {
        return get_normality_statistic();
    }
    return get_aux_normality_statistic(rff_feature_levels.size() - 1);
}

double mvn_test::get_aux_normality_statistic(size_t level) const {
    if (!has_aux_statistic()) {
        return get_normality_statistic();
    }
    if (level >= rff_feature_levels.size() || level >= rff_pairwise_stats.size()) {
        throw std::logic_error("RFF ladder state is inconsistent");
    }

    if (effect_size <= dimensions()) {
        throw std::logic_error("Too few points.");
    }
    const double beta = betas.front();
    double stat = std::pow(1 + 2 * std::pow(beta, 2), dimensions() / -2.0);
    stat += (1.0 / ((double)effect_size * effect_size)) * rff_pairwise_stats[level];
    stat -= (2.0 / (effect_size * std::pow(1 + std::pow(beta, 2.0), dimensions() / 2.0))) * rff_center_stat;
    return stat;
}

double mvn_test::aux_delta_last_swap(size_t level) const {
    if (!has_aux_statistic()) {
        throw std::logic_error("Auxiliary delta requested but hybrid mode is disabled");
    }
    if (!last_swap_has_equal_effect_size()) {
        throw std::logic_error("Auxiliary delta requires unchanged effect size (equal cluster sizes)");
    }
    if (level >= rff_feature_levels.size()) {
        throw std::logic_error("RFF ladder level is out of range");
    }
    if (effect_size <= dimensions()) {
        throw std::logic_error("Too few points.");
    }

    const size_t a = (size_t)latest_subset_point;
    const size_t b = (size_t)latest_replacing_point;
    const size_t dim = rff_feature_levels[level];
    const float* phi_a = rff_stats->features_ptr(a);
    const float* phi_b = rff_stats->features_ptr(b);
    const double beta = betas.front();
    const double n = (double)effect_size;
    const double denom_center = std::pow(1 + std::pow(beta, 2.0), dimensions() / 2.0);

    const double dot_sum_b = dot_float_sse(rff_subset_feature_sum.data(), phi_b, dim);
    const double dot_sum_a = dot_float_sse(rff_subset_feature_sum.data(), phi_a, dim);
    const double self_b = dot_float_sse(phi_b, phi_b, dim);
    const double self_a = dot_float_sse(phi_a, phi_a, dim);
    const double cross_ab = dot_float_sse(phi_a, phi_b, dim);

    const double delta_pairwise = 2.0 * (dot_sum_b - dot_sum_a) - (self_b + self_a - 2.0 * cross_ab);
    const double delta_center = rff_stats->centered_stat(b) - rff_stats->centered_stat(a);

    double delta_stat = (1.0 / (n * n)) * delta_pairwise;
    delta_stat -= (2.0 / (n * denom_center)) * delta_center;
    return delta_stat;
}

bool mvn_test::last_swap_has_equal_effect_size() const {
    if (latest_subset_point < 0 || latest_replacing_point < 0) {
        return false;
    }
    return clustering->cluster_size((size_t)latest_subset_point) == clustering->cluster_size((size_t)latest_replacing_point);
}

namespace {

double exact_centered_cluster(const mvn::mahalanobis_distances& distances, const mvn::Clustering& clst, size_t cluster, double k_center) {
    double acc = 0.0;
    for (int el : clst.elements(cluster)) {
        acc += std::exp(k_center * distances.distance((unsigned)el));
    }
    return acc;
}

double exact_pairwise_clusters(const mvn::mahalanobis_distances& distances, const mvn::Clustering& clst, size_t a, size_t b, double k_pw) {
    double acc = 0.0;
    const auto& lhs = clst.elements(a);
    const auto& rhs = clst.elements(b);
    for (int i : lhs) {
        for (int j : rhs) {
            acc += std::exp(k_pw * distances.interpoint_distance((unsigned)i, (unsigned)j));
        }
    }
    if (a == b) {
        acc *= 0.5;
    }
    return acc;
}

}

double mvn_test::exact_delta_last_swap() const {
    if (!last_swap_has_equal_effect_size()) {
        throw std::logic_error("Exact delta requires unchanged effect size (equal cluster sizes)");
    }
    if (betas.empty()) {
        throw std::logic_error("No beta configured");
    }
    const double beta = betas.front();
    const double k_center = -(beta * beta / (2 * (1 + beta * beta)));
    const double k_pw = -(beta * beta) / 2.0;

    const size_t a = (size_t)latest_subset_point;
    const size_t b = (size_t)latest_replacing_point;
    const double n = (double)effect_size;
    const double denom_center = std::pow(1 + std::pow(beta, 2.0), dimensions() / 2.0);

    double sum_bS = 0.0;
    double sum_aS = 0.0;
    for (size_t s : subset) {
        if (s == b) {
            continue;
        }
        sum_bS += exact_pairwise_clusters(*distances, *clustering, b, s, k_pw);
        sum_aS += exact_pairwise_clusters(*distances, *clustering, a, s, k_pw);
    }
    const double self_b = exact_pairwise_clusters(*distances, *clustering, b, b, k_pw);
    const double self_a = exact_pairwise_clusters(*distances, *clustering, a, a, k_pw);

    const double delta_pairwise = 2.0 * (sum_bS - sum_aS) + (self_b - self_a);
    const double delta_center = exact_centered_cluster(*distances, *clustering, b, k_center) -
                                exact_centered_cluster(*distances, *clustering, a, k_center);

    double delta_stat = (1.0 / (n * n)) * delta_pairwise;
    delta_stat -= (2.0 / (n * denom_center)) * delta_center;
    return delta_stat;
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

double mvn_test::get_normality_statistic() {
    return static_cast<const mvn_test&>(*this).get_normality_statistic();
}

const std::vector<size_t>& mvn_test::current_subset() const {
    return subset;
}

mvn_stats::mvn_stats(const mahalanobis_distances& distances, const Clustering& clst, double beta,
                     bool use_nystrom, size_t n_features,
                     bool use_rff, size_t rff_features, uint32_t seed)
    :mahalanobis_centered(clst.size()),
     mahalanobis_pairwise() {
    this->distances = &distances;
    this->clustering = &clst;
    const size_t n_clusters = clst.size();

    double k_center = -(beta * beta / (2 * (1 + beta * beta)));
    k_pw = -(beta * beta) / 2.0;

    for (size_t cl = 0; cl < n_clusters; cl++) {
        for (int el: clst.elements(cl)) {
            mahalanobis_centered[cl] += std::exp(k_center * distances.distance(el));
        }
    }

    size_t n_samples = 0;
    for (size_t cl = 0; cl < n_clusters; ++cl) {
        n_samples += clst.elements(cl).size();
    }

    size_t feature_count = 0;
    if (use_rff) {
        feature_mode = true;
        feature_dim = std::max<size_t>(1, rff_features);
    } else {
        const size_t requested_features = (n_features == 0) ? NYSTROM_DEFAULT_MAX_FEATURES : n_features;
        feature_count = std::min(n_samples, requested_features);
        feature_mode = use_nystrom && (feature_count > 0);
        feature_dim = feature_count;
    }

    if (use_rff) {
        std::vector<int> sample_cluster_map(n_samples, -1);
        for (size_t cl = 0; cl < n_clusters; ++cl) {
            for (int el : clst.elements(cl)) {
                sample_cluster_map[el] = static_cast<int>(cl);
            }
        }

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(distances.inv_cov());
        if (eig.info() != Eigen::Success) {
            throw std::logic_error("Failed to decompose inv_cov for RFF");
        }
        Eigen::VectorXd evals = eig.eigenvalues();
        Eigen::MatrixXd evecs = eig.eigenvectors();
        for (int i = 0; i < evals.size(); ++i) {
            evals(i) = std::max(evals(i), 0.0);
        }
        Eigen::MatrixXd sqrt_inv = evecs * evals.array().sqrt().matrix().asDiagonal() * evecs.transpose();

        std::mt19937 rng(seed);
        std::normal_distribution<double> normal(0.0, 1.0);
        constexpr double PI = 3.141592653589793238462643383279502884;
        std::uniform_real_distribution<double> unif(0.0, 2.0 * PI);

        Eigen::MatrixXd W(distances.data().rows(), (Eigen::Index)feature_dim);
        for (Eigen::Index d = 0; d < (Eigen::Index)feature_dim; ++d) {
            Eigen::VectorXd z(distances.data().rows());
            for (Eigen::Index r = 0; r < z.size(); ++r) {
                z(r) = normal(rng);
            }
            Eigen::VectorXd w = beta * (sqrt_inv * z);
            W.col(d) = w;
        }
        std::vector<double> b(feature_dim);
        for (size_t d = 0; d < feature_dim; ++d) {
            b[d] = unif(rng);
        }
        const double scale = std::sqrt(2.0 / (double)feature_dim);

        cluster_features.assign(n_clusters, std::vector<float>(feature_dim, 0.0f));
        for (size_t sample = 0; sample < n_samples; ++sample) {
            int cl = sample_cluster_map[sample];
            if (cl < 0) {
                continue;
            }
            Eigen::VectorXd x = distances.data().col((Eigen::Index)sample);
            Eigen::RowVectorXd proj = x.transpose() * W;
            auto& dst = cluster_features[(size_t)cl];
            for (size_t d = 0; d < feature_dim; ++d) {
                dst[d] += static_cast<float>(scale * std::cos(proj((Eigen::Index)d) + b[d]));
            }
        }

        mahalanobis_pairwise.clear();
        return;
    }

    if (feature_mode) {
        std::vector<int> sample_cluster_map(n_samples, -1);
        for (size_t cl = 0; cl < n_clusters; ++cl) {
            for (int el : clst.elements(cl)) {
                sample_cluster_map[el] = static_cast<int>(cl);
            }
        }

        std::mt19937 local_wheel(seed);
        std::vector<size_t> landmarks = choose_landmarks(n_samples, feature_count, local_wheel);

        Eigen::MatrixXd W(feature_count, feature_count);
        for (size_t i = 0; i < feature_count; ++i) {
            for (size_t j = i; j < feature_count; ++j) {
                double val = std::exp(k_pw * distances.interpoint_distance(landmarks[i], landmarks[j]));
                W(i, j) = val;
                W(j, i) = val;
            }
        }

        W = 0.5 * (W + W.transpose());
        const double mean_diag = std::max(std::abs(W.diagonal().mean()), 1.0);
        const double ridge = std::max(NYSTROM_RIDGE_MIN, NYSTROM_RIDGE_SCALE * mean_diag);
        W.diagonal().array() += ridge;

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(W);
        if (eig.info() != Eigen::Success) {
            feature_mode = false;
        } else {
            Eigen::VectorXd evals = eig.eigenvalues();
            Eigen::MatrixXd evecs = eig.eigenvectors();

            const double max_eval = std::max(evals.maxCoeff(), ridge);
            const double floor_rel = NYSTROM_REL_EIGEN_FLOOR * max_eval;
            const double floor_cond = max_eval / NYSTROM_MAX_COND;
            const double eigen_floor = std::max({ridge, floor_rel, floor_cond});
            for (int i = 0; i < evals.size(); ++i) {
                evals(i) = std::max(evals(i), eigen_floor);
            }

            Eigen::VectorXd inv_sqrt = evals.array().sqrt().inverse();
            Eigen::MatrixXd transform = evecs * inv_sqrt.asDiagonal();

            cluster_features.assign(n_clusters, std::vector<float>(feature_count, 0.0f));
            Eigen::RowVectorXd c_row(feature_count);
            for (size_t sample = 0; sample < n_samples; ++sample) {
                for (size_t l = 0; l < feature_count; ++l) {
                    c_row(static_cast<Eigen::Index>(l)) = std::exp(k_pw * distances.interpoint_distance(sample, landmarks[l]));
                }
                Eigen::RowVectorXd feature = c_row * transform;
                int cl = sample_cluster_map[sample];
                if (cl < 0) {
                    continue;
                }
                auto& dst = cluster_features[static_cast<size_t>(cl)];
                for (size_t d = 0; d < feature_count; ++d) {
                    dst[d] += static_cast<float>(feature(static_cast<Eigen::Index>(d)));
                }
            }
        }
    }

    if (!feature_mode) {
        cluster_features.clear();
        feature_dim = 0;
        mahalanobis_pairwise.resize(n_clusters);
        for (size_t cl = 0; cl < n_clusters; ++cl) {
            mahalanobis_pairwise[cl].resize(n_clusters);
            for (size_t pair_cl = 0; pair_cl < n_clusters; ++pair_cl) {
                double acc = 0.0;
                for (int el : clst.elements(cl)) {
                    for (int pair_el : clst.elements(pair_cl)) {
                        const double value = std::exp(k_pw * distances.interpoint_distance((unsigned)el, (unsigned)pair_el));
                        acc += (cl == pair_cl) ? 0.5 * value : value;
                    }
                }
                mahalanobis_pairwise[cl][pair_cl] = acc;
            }
        }
        return;
    }

    mahalanobis_pairwise.resize(n_clusters);
    for (size_t cl = 0; cl < n_clusters; ++cl) {
        mahalanobis_pairwise[cl].resize(n_clusters);
        for (size_t pair_cl = 0; pair_cl < n_clusters; ++pair_cl) {
            mahalanobis_pairwise[cl][pair_cl] = feature_pairwise_stat(cl, pair_cl);
        }
    }
}

double mvn_stats::feature_pairwise_stat(size_t i, size_t j) const {
    const auto& lhs = cluster_features[i];
    const auto& rhs = cluster_features[j];
    double dot = dot_float_sse(lhs.data(), rhs.data(), feature_dim);
    if (i == j) {
        return 0.5 * dot;
    }
    return dot;
}

double mvn_stats::pairwise_stat(size_t i, size_t j) const {
    if (feature_mode) {
        return feature_pairwise_stat(i, j);
    }

    if (!mahalanobis_pairwise.empty()) {
        return mahalanobis_pairwise[i][j];
    }

    if (!distances || !clustering) {
        throw std::logic_error("Exact pairwise requested but distances/clustering are not set");
    }

    double acc = 0.0;
    const auto& lhs = clustering->elements(i);
    const auto& rhs = clustering->elements(j);
    for (int a : lhs) {
        for (int b : rhs) {
            acc += std::exp(k_pw * distances->interpoint_distance((unsigned)a, (unsigned)b));
        }
    }
    if (i == j) {
        acc *= 0.5;
    }
    return acc;
}

double mvn_stats::centered_stat(size_t i) const {
    return mahalanobis_centered[i];
}

double mvn_stats::sum_pairwise(size_t point, const std::vector<size_t>& ss) const {
    if (feature_mode) {
        double ret = 0.0;
        for (auto s: ss) {
            ret += feature_pairwise_stat(point, s);
        }
        return ret;
    }

    if (!mahalanobis_pairwise.empty()) {
        double ret = 0.0;
        const auto& pairwise_row = mahalanobis_pairwise[point];
        for (auto s: ss) {
            ret += pairwise_row[s];
        }
        return ret;
    }

    double ret = 0.0;
    for (auto s: ss) {
        ret += pairwise_stat(point, s);
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
        if (latest_subset_point == -1) {
            throw std::logic_error("Unable to reject the last swap.");
        }
        replacing_point = latest_subset_point;
        latest_subset_point = -1;
        latest_replacing_point = -1;
    } else {
        std::uniform_int_distribution<unsigned> subset_unif(0, subsample_size() - 1);
        std::swap(subset[subset_unif(wheel)], subset.back());
        replacing_point = sampler.sample();
    }

    auto subset_point = subset.back();
    if (!reject_last) {
        latest_subset_point = subset_point;
        latest_replacing_point = replacing_point;
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
    subset_feature_sum(other.subset_feature_sum),
     betas(other.betas),
     clustering(other.clustering),
     p(other.p),
     n(other.n),
     effect_size(other.effect_size),
     latest_subset_point(other.latest_subset_point),
    latest_replacing_point(other.latest_replacing_point),
    use_nystrom(other.use_nystrom),
    n_features(other.n_features),
    use_hybrid(other.use_hybrid),
    rff_stats(other.rff_stats),
    rff_feature_levels(other.rff_feature_levels),
    rff_pairwise_stats(other.rff_pairwise_stats),
    rff_center_stat(other.rff_center_stat),
    rff_subset_feature_sum(other.rff_subset_feature_sum),
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
    return clusters[i];
}

size_t Clustering::cluster_size(size_t i) const {
    return clusters[i].size();
}

mahalanobis_distances::mahalanobis_distances(std::shared_ptr<const Matrix> X, const Matrix& cov, const Vector& mean)
        :X(std::move(X)) {
    if (!this->X || this->X->cols() == 0 || this->X->rows() == 0) {
        throw std::invalid_argument("Matrix is empty");
    }

    Eigen::FullPivHouseholderQR<Matrix> qr(cov);
    if (!qr.isInvertible()) {
        throw std::logic_error("Non-invertible matrix. Must not happen.");
    }
    S_inv = qr.inverse();

    S_inv_X = S_inv * (*this->X);
    S_inv_mean = S_inv * mean;
    mu_mu = mean.dot(S_inv_mean);

    const Eigen::Index n = this->X->cols();
    quad_x.assign((size_t)n, 0.0);
    x_mu.assign((size_t)n, 0.0);
    for (Eigen::Index i = 0; i < n; ++i) {
        const auto xi = this->X->col(i);
        const auto yi = S_inv_X.col(i);
        quad_x[(size_t)i] = xi.dot(yi);
        x_mu[(size_t)i] = xi.dot(S_inv_mean);
    }
}

double mahalanobis_distances::interpoint_distance(unsigned i, unsigned j) const {
    const double xixj = (*X).col((Eigen::Index)i).dot(S_inv_X.col((Eigen::Index)j));
    return quad_x[i] - 2.0 * xixj + quad_x[j];
}

double mahalanobis_distances::distance(unsigned i) const {
    return quad_x[i] - 2.0 * x_mu[i] + mu_mu;
}

void mvn_test::remove(unsigned point) {
    for (size_t i = 0; i < stats.size(); i++) {
        if (stats[i]->is_feature_mode()) {
            const float* phi = stats[i]->features_ptr(point);
            auto& sum = subset_feature_sum[i];
            const size_t dim = stats[i]->features_dim();
            const double dot = dot_float_sse(phi, sum.data(), dim);
            const double self = dot_float_sse(phi, phi, dim);
            const double sp = dot - 0.5 * self;
            pairwise_stat[i] -= 2.0 * sp;
            for (size_t d = 0; d < dim; ++d) {
                sum[d] -= phi[d];
            }
        } else {
            pairwise_stat[i] -= 2.0 * stats[i]->sum_pairwise(point, subset);
        }
    }

    for (size_t i = 0; i < stats.size(); i++) {
        center_stat[i] -= stats[i]->centered_stat(point);
    }

    if (has_aux_statistic()) {
        const float* phi = rff_stats->features_ptr(point);
        for (size_t level = 0; level < rff_feature_levels.size(); ++level) {
            const size_t dim = rff_feature_levels[level];
            const double dot = dot_float_sse(phi, rff_subset_feature_sum.data(), dim);
            const double self = dot_float_sse(phi, phi, dim);
            const double sp = dot - 0.5 * self;
            rff_pairwise_stats[level] -= 2.0 * sp;
        }
        const size_t full_dim = rff_stats->features_dim();
        for (size_t d = 0; d < full_dim; ++d) {
            rff_subset_feature_sum[d] -= phi[d];
        }
        rff_center_stat -= rff_stats->centered_stat(point);
    }
}

void mvn_test::add(unsigned int point) {
    for (size_t i = 0; i < stats.size(); i++) {
        if (stats[i]->is_feature_mode()) {
            const float* phi = stats[i]->features_ptr(point);
            auto& sum = subset_feature_sum[i];
            const size_t dim = stats[i]->features_dim();
            for (size_t d = 0; d < dim; ++d) {
                sum[d] += phi[d];
            }
            const double dot = dot_float_sse(phi, sum.data(), dim);
            const double self = dot_float_sse(phi, phi, dim);
            const double sp = dot - 0.5 * self;
            pairwise_stat[i] += 2.0 * sp;
        } else {
            pairwise_stat[i] += 2.0 * stats[i]->sum_pairwise(point, subset);
        }
    }
    for (size_t i = 0; i < stats.size(); i++) {
        center_stat[i] += stats[i]->centered_stat(point);
    }

    if (has_aux_statistic()) {
        const float* phi = rff_stats->features_ptr(point);
        const size_t full_dim = rff_stats->features_dim();
        for (size_t d = 0; d < full_dim; ++d) {
            rff_subset_feature_sum[d] += phi[d];
        }
        for (size_t level = 0; level < rff_feature_levels.size(); ++level) {
            const size_t dim = rff_feature_levels[level];
            const double dot = dot_float_sse(phi, rff_subset_feature_sum.data(), dim);
            const double self = dot_float_sse(phi, phi, dim);
            const double sp = dot - 0.5 * self;
            rff_pairwise_stats[level] += 2.0 * sp;
        }
        rff_center_stat += rff_stats->centered_stat(point);
    }
}

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
    return active_tree[el_pos(n)];
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
    active_tree[pos] = true;
    update(pos);
}

size_t RandomSampler::sample() {
    if (size < 2) {
        throw std::logic_error("Too little segment tree");
    }
    size_t node = 0;
    while (!is_leaf(node)) {
        auto chld = children(node);
        if (active_tree[chld.first] == 0) {
            if (active_tree[chld.second] == 0) {
                throw std::logic_error("No active elements in subtree");
            }
            node = chld.second;
        } else if(active_tree[chld.second] == 0) {
            node = chld.first;
        } else {
            double l = segment_tree[chld.first];
            double r = segment_tree[chld.second];
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
    if (!active_tree[node]) {
        throw std::logic_error("Sampled element is not active.");
    }
    int element = node - aug_nodes;
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
    active_tree[node] = active_tree[chs.first] + active_tree[chs.second];
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
    return active_tree[0];
}

}
