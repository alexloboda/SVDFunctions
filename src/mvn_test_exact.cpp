#include <cmath>
#include <cctype>
#include <limits>
#include <stdexcept>

#include "include/third-party/cxxpool.h"

#include "include/mvn_test.h"

namespace mvn {

namespace {

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

std::string to_lower_ascii(std::string s) {
    for (char& ch : s) {
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    }
    return s;
}

} // namespace

mvn_test_method parse_mvn_test_method(const std::string& method) {
    const std::string m = to_lower_ascii(method);
    if (m == "exact" || m == "kernel" || m == "full") {
        return mvn_test_method::exact;
    }
    if (m == "rff" || m == "approx" || m == "approximate") {
        return mvn_test_method::rff;
    }
    throw std::invalid_argument("Unknown mvn test method: " + method);
}

std::shared_ptr<mvn_test_base> make_mvn_test(std::shared_ptr<const Matrix> X,
                                            const Clustering& clst,
                                            const Vector& mean,
                                            const Matrix& cov,
                                            mvn_test_method method,
                                            size_t rff_dim) {
    switch (method) {
        case mvn_test_method::exact:
            return std::make_shared<mvn_test_exact>(std::move(X), clst, cov, mean);
        case mvn_test_method::rff:
            return std::make_shared<mvn_test_rff>(std::move(X), clst, cov, mean, rff_dim);
        default:
            throw std::invalid_argument("Unsupported mvn test method");
    }
}

class mvn_test_exact::mahalanobis_distances {
    std::vector<std::vector<double>> inter;
    std::vector<double> dist;

public:
    mahalanobis_distances(std::shared_ptr<const Matrix> X, const Matrix& cov, const Vector& mean)
            :dist(static_cast<size_t>(X->cols())) {
        if (X->cols() == 0 || X->rows() == 0) {
            throw std::invalid_argument("Matrix is empty");
        }
        if (cov.rows() != X->rows() || cov.cols() != X->rows()) {
            throw std::invalid_argument("Covariance dimension mismatch");
        }
        if (mean.size() != X->rows()) {
            throw std::invalid_argument("Mean vector dimension mismatch");
        }

        inter.resize(static_cast<size_t>(X->cols()));

        Eigen::FullPivHouseholderQR<Matrix> qr(cov);
        if (!qr.isInvertible()) {
            throw std::logic_error("Non-invertible matrix. Must not happen.");
        }
        Matrix cov_inv = qr.inverse();

        Vector ximu = X->transpose() * cov_inv * mean;
        Vector muxi = mean.transpose() * cov_inv * *X;
        Matrix distances = X->transpose() * cov_inv * *X;
        double mumu = mean.transpose() * cov_inv * mean;

        std::vector<double> diag(static_cast<size_t>(distances.rows()));
        for (Eigen::Index i = 0; i < distances.rows(); i++) {
            diag[static_cast<size_t>(i)] = distances(i, i);
        }

        for (Eigen::Index i = 0; i < X->cols(); i++) {
            inter[static_cast<size_t>(i)].resize(static_cast<size_t>(X->cols()));
            dist[static_cast<size_t>(i)] = diag[static_cast<size_t>(i)] - ximu(i) - muxi(i) + mumu;
            for (Eigen::Index j = 0; j < X->cols(); j++) {
                inter[static_cast<size_t>(i)][static_cast<size_t>(j)] =
                        diag[static_cast<size_t>(i)] - 2 * distances(i, j) + diag[static_cast<size_t>(j)];
            }
        }
    }

    double distance(unsigned i) const {
        return dist.at(i);
    }

    double interpoint_distance(unsigned i, unsigned j) const {
        return inter.at(i).at(j);
    }
};

class mvn_test_exact::mvn_stats {
    std::vector<double> mahalanobis_centered;
    std::vector<std::vector<double>> mahalanobis_pairwise;

public:
    mvn_stats(const mahalanobis_distances& distances, const Clustering& clst, double beta)
        :mahalanobis_centered(clst.size()),
         mahalanobis_pairwise() {
        mahalanobis_pairwise.resize(clst.size());
        size_t n = clst.size();

        double k_center = -(beta * beta / (2 * (1 + beta * beta)));
        double k_pw = -(beta * beta) / 2.0;

        for (size_t cl = 0; cl < n; cl++) {
            mahalanobis_pairwise[cl].resize(clst.size());
            for (int el: clst.elements(cl)) {
                mahalanobis_centered[cl] += std::exp(k_center * distances.distance(static_cast<unsigned>(el)));
                for (size_t pair_cl = 0; pair_cl < n; pair_cl++) {
                    for (int pair_el: clst.elements(pair_cl)) {
                        const double v = std::exp(k_pw * distances.interpoint_distance(static_cast<unsigned>(el),
                                                                                      static_cast<unsigned>(pair_el)));
                        if (cl == pair_cl) {
                            mahalanobis_pairwise[cl][pair_cl] += 0.5 * v;
                        } else {
                            mahalanobis_pairwise[cl][pair_cl] += v;
                        }
                    }
                }
            }
        }
    }

    double sum_pairwise(size_t point, const std::vector<size_t>& ss) const {
        double ret = 0.0;
        for (auto s: ss) {
            ret += mahalanobis_pairwise.at(point).at(s);
        }
        return ret;
    }

    double centered_stat(size_t i) const {
        return mahalanobis_centered.at(i);
    }
};

mvn_test_exact::mvn_test_exact(std::shared_ptr<const Matrix> X,
                               const Clustering& clst,
                               const Matrix& cov,
                               const Vector& mean)
        :distances{std::make_shared<mahalanobis_distances>(X, cov, mean)},
         clustering(std::make_shared<Clustering>(clst)),
         wheel(std::random_device()()) {
    if (X->cols() == 0 || X->rows() == 0) {
        throw std::invalid_argument("Matrix is empty");
    }

    betas = {0.8};

    pairwise_stat.resize(betas.size(), 0.0);
    center_stat.resize(betas.size(), 0.0);

    n = clst.size();
    p = static_cast<size_t>(X->rows());
    effect_size = 0;
    latest_subset_point = -1;

    if (n == 0) {
        throw std::invalid_argument("Clustering is empty");
    }
    if (p == 0) {
        throw std::invalid_argument("Matrix has zero rows");
    }
    if (n <= p) {
        throw std::logic_error("Too few points.");
    }

    cxxpool::thread_pool pool(std::thread::hardware_concurrency());
    std::vector<std::future<std::shared_ptr<mvn_stats>>> futures;
    futures.reserve(betas.size());
    for (double beta: betas) {
        futures.push_back(pool.push([this, clst, beta]() -> std::shared_ptr<mvn_stats> {
            return std::make_shared<mvn_stats>(*distances, clst, beta);
        }));
    }
    for (size_t i = 0; i < betas.size(); i++) {
        stats.push_back(futures[i].get());
    }

    std::vector<double> lls;
    lls.reserve(n);
    for (size_t i = 0; i < n; i++) {
        auto ids = clst.elements(i);
        (void)loglikelihood(ids);
        lls.push_back(0.0);
    }
    sampler = RandomSampler(lls, wheel());

    while (effect_size < p + 1) {
        add_one();
    }
}

mvn_test_exact::mvn_test_exact(const mvn_test_exact& other)
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
    // NOTE: don't read other's RNG state here; copy/clone may run in parallel.
    wheel{std::random_device()()},
     subset(other.subset) {}

size_t mvn_test_exact::dimensions() const {
    return p;
}

size_t mvn_test_exact::subsample_size() const {
    return subset.size();
}

size_t mvn_test_exact::sample_size() const {
    return n;
}

const std::vector<size_t>& mvn_test_exact::current_subset() const {
    return subset;
}

double mvn_test_exact::get_normality_statistic() {
    if (effect_size <= dimensions()) {
        throw std::logic_error("Too few points.");
    }
    double max_stat = 0.0;
    for (size_t i = 0; i < betas.size(); i++) {
        double stat = std::pow(1 + 2 * std::pow(betas[i], 2), dimensions() / -2.0);
        stat += (1.0 / (static_cast<double>(effect_size) * effect_size)) * pairwise_stat[i];
        stat -= (2.0 / (effect_size * std::pow(1 + std::pow(betas[i], 2.0), dimensions() / 2.0))) * center_stat[i];
        max_stat = std::max(max_stat, stat);
    }
    return max_stat;
}

void mvn_test_exact::swap_once(bool reject_last) {
    if (subset.empty() || sampler.n_active() == 0) {
        throw std::logic_error("Unable to swap points.");
    }

    int replacing_point = -1;
    if (reject_last) {
        if (latest_subset_point == -1) {
            throw std::logic_error("No last point to reject");
        }
        replacing_point = latest_subset_point;
        latest_subset_point = -1;
    } else {
        std::uniform_int_distribution<unsigned> subset_unif(0, subsample_size() - 1);
        std::swap(subset[subset_unif(wheel)], subset.back());
        replacing_point = static_cast<int>(sampler.sample());
    }

    auto subset_point = subset.back();
    if (!reject_last) {
        latest_subset_point = static_cast<int>(subset_point);
    }

    effect_size += clustering->cluster_size(static_cast<size_t>(replacing_point)) - clustering->cluster_size(subset_point);

    remove(static_cast<unsigned>(subset_point));

    sampler.enable(subset.back());
    sampler.disable(static_cast<size_t>(replacing_point));

    subset.pop_back();
    subset.push_back(static_cast<size_t>(replacing_point));

    add(static_cast<unsigned>(replacing_point));
}

void mvn_test_exact::add_one() {
    if (sampler.n_active() == 0) {
        throw std::logic_error("Can't add point to the model.");
    }
    latest_subset_point = -1;

    size_t point = sampler.sample();
    effect_size += clustering->cluster_size(point);
    sampler.disable(point);

    subset.push_back(point);
    add(static_cast<unsigned>(point));
}

void mvn_test_exact::remove(unsigned point) {
    for (size_t i = 0; i < stats.size(); i++) {
        pairwise_stat[i] -= 2 * stats[i]->sum_pairwise(point, subset);
    }
    for (size_t i = 0; i < stats.size(); i++) {
        center_stat[i] -= stats[i]->centered_stat(point);
    }
}

void mvn_test_exact::add(unsigned point) {
    for (size_t i = 0; i < stats.size(); i++) {
        pairwise_stat[i] += 2 * stats[i]->sum_pairwise(point, subset);
    }
    for (size_t i = 0; i < stats.size(); i++) {
        center_stat[i] += stats[i]->centered_stat(point);
    }
}

std::shared_ptr<mvn_test_base> mvn_test_exact::clone() const {
    return std::make_shared<mvn_test_exact>(*this);
}

std::vector<double> mvn_test_exact::loglikelihood(const std::vector<int>& ids) const {
    std::vector<double> ret;
    ret.reserve(ids.size());
    for (int id: ids) {
        double dist = std::sqrt(distances->distance(static_cast<unsigned>(id)));
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
