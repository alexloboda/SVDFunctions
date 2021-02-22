#include "include/matching.h"
#include "include/hw.h"
#include "include/lm.h"

#include <numeric>
#include <utility>

namespace {

double lambda_1000(double lambda_observed, int cases, int controls) {
    double coefficient = 1.0 / cases + 1.0 / controls;
    return 1.0 + ((lambda_observed - 1) * coefficient) / (2 * 0.001);
}

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

}

namespace matching {

matching::matching(Eigen::MatrixXi controls_gmatrix,
                   Eigen::MatrixXd controls_space,
                   Clustering clustering) : controls_gmatrix(std::move(controls_gmatrix)),
                                            controls_space(std::move(controls_space)),
                                            clustering(std::move(clustering)) {}

void matching::set_soft_threshold(lambda_range range) {
    soft_threshold = range;
}

void matching::set_hard_threshold(lambda_range range) {
    hard_threshold = range;
}

std::vector<lm> init_lms(const std::vector<Counts>& case_counts) {
    std::vector<Counts> counts(case_counts.size());
    std::vector<lm> lms;
    for (size_t i = 0; i < case_counts.size(); i++) {
        lm model(6);
        for (size_t j = 0; j < 3; j++) {
            model.set(3 + (int)j, j, 1, case_counts[i][j]);
        }
        lms.push_back(model);
    }

    return lms;
}

matching_results matching::match(const std::vector<Counts>& case_counts, unsigned min_controls) {
    size_t n_variants = case_counts.size();

    std::vector<bool> snp_mask = check_user_counts(case_counts);
    auto lms = init_lms(case_counts);

    double lambda = std::numeric_limits<double>::infinity();

    std::vector<double> lambdas;

    std::vector<int> optimal_controls;
    std::vector<int> lambda_i;
    std::vector<int> pvals_num;
    std::vector<double> optimal_pvals;

    for (size_t k = subsampling.min_size(), step = 0; k <= subsampling.max_size(); k++, step++) {
        if (step % 100 == 0) {
            interrupts_checker();
        }

        std::vector<double> pvals;

        auto control_groups = subsampling.get_solution(k);
        std::vector<int> controls;
        for (int group: control_groups) {
            auto group_elements = clustering.elements(group);
            controls.insert(controls.end(), group_elements.begin(), group_elements.end());
        }

        for (size_t j = 0; j < n_variants; j++) {
            if (!snp_mask[j]) {
                continue;
            }

            Counts controls_counts = count_controls(controls, j);

            if (!check_counts(controls_counts[0], controls_counts[1], controls_counts[2])) {
                continue;
            }

            for (int i = 0; i < 3; i++) {
                lms[j].set(i, i, 0, controls_counts[i]);
            }

            lms[j].solve();
            int rank = controls_counts.sum() + case_counts[j].sum();
            pvals.push_back(pval_t(lms[j].compute_t(rank - 2), rank - 2));
        }

        if (controls.size() >= min_controls && pvals.size() > 10) {
            double cur_lambda = get_lambda(pvals);

            int cases = std::max_element(case_counts.begin(), case_counts.end(), [](const Counts& a, const Counts& b) {
                return a.sum() < b.sum();
            })->sum();
            auto num_of_controls = controls.size();
            cur_lambda = lambda_1000(cur_lambda, cases, num_of_controls);

            lambdas.push_back(cur_lambda);
            lambda_i.push_back(controls.size());
            pvals_num.push_back(pvals.size());

            if (hard_threshold.in(cur_lambda)) {
                if (soft_threshold.in(cur_lambda) || soft_threshold.distance(cur_lambda) < soft_threshold.distance(lambda)) {
                    optimal_controls = controls;
                    lambda = cur_lambda;
                    optimal_pvals = pvals;
                }
            }
        }
    }
    return {std::move(optimal_controls), std::move(optimal_pvals), std::move(lambdas), std::move(lambda_i), std::move(pvals_num)};
}

void matching::process_mvn(const Matrix& cov, const Vector& mean) {
    Eigen::MatrixXd controls_gmatrix_float(controls_space.rows(), clustering.size());
    for (size_t i = 0; i < clustering.size(); i++) {
        controls_gmatrix_float.col(i) = cluster_mean(i);
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(controls_space);
    const Matrix& U = svd.matrixU();

    Vector rs_mean = U.transpose() * mean;
    Matrix rs_cov = U.transpose() * cov;

    Matrix rs_controls_gmatrix = U.transpose() * controls_gmatrix_float;
    subsampling = mvn::subsample(rs_controls_gmatrix, rs_mean, rs_cov);
    subsampling.run(100000, 10, 1.0 / 3.0);
}

Eigen::VectorXd matching::cluster_mean(int cluster) const {
    auto objects = clustering.elements(cluster);
    return std::accumulate(objects.begin(), objects.end(), Vector(controls_space.rows()), [&](Vector sum, int obj) {
        return sum + controls_space.col(obj);
    });
}

Counts matching::count_controls(const std::vector<int>& controls, size_t variant) {
    Counts counts;
    for (int sample: controls) {
        double value = controls_gmatrix(variant, sample);
        if (!std::isnan(value)) {
            ++counts[controls_gmatrix(variant, sample)];
        }
    }
    return counts;
}

void matching::set_qchi_sq_function(const std::function<double(double)>& f) {
    qchisq = f;
}

double matching::get_lambda(std::vector<double>& pvals) {
    int n_pvals = pvals.size();

    lm pvals_lm(n_pvals);
    std::sort(pvals.begin(), pvals.end());
    for (size_t j = 0; j < pvals.size(); j++) {
        pvals_lm.set(j, chi2(j, n_pvals, qchisq), qchisq(pvals[j]), 1, false);
    }
    pvals_lm.solve();

    return pvals_lm.get_lambda();
}

lambda_range::lambda_range(double lb, double ub) : lb(lb), ub(ub) {}

lambda_range::lambda_range() : lb(0.0), ub(std::numeric_limits<double>::infinity()) {}

bool lambda_range::in(double lambda) {
    return lambda >= lb && lambda < ub;
}

double lambda_range::distance(double lambda) {
    if (in(lambda)) {
        throw std::invalid_argument("Lambda is in range");
    }

    return std::max(lambda - ub, lb - lambda);
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
    return n;
}

std::vector<int> Clustering::elements(size_t i) const {
    return clusters.at(i);
}

matching_results::matching_results(std::vector<int>&& prefix, std::vector<double>&& p_values,
                                   std::vector<double>&& lmbds, std::vector<int>&& lmbd_i,
                                   std::vector<int>&& pvals_number)
        : optimal_prefix(std::move(prefix)), pvals(std::move(p_values)),
          lambdas(std::move(lmbds)), lambda_i(std::move(lmbd_i)), pvals_num(std::move(pvals_number)) {}

}