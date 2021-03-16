#include "include/matching.h"
#include "include/hw.h"
#include "include/lm.h"

#include <numeric>
#include <utility>
#include <fstream>
#include <Rcpp.h>

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
                   mvn::Clustering clustering) : controls_gmatrix(std::move(controls_gmatrix)),
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

void matching::process_mvn(const Matrix& directions, Vector mean, int threads) {
    Rcpp::Rcerr << "Starting processing controls space." << std::endl;
    Rcpp::Rcerr << "The size of controls space is " << controls_space.rows() << " by " << controls_space.cols() << std::endl;

    Vector controls_mean = controls_space.rowwise().mean();
    Rcpp::Rcerr << "Controls mean size: " << controls_mean.size() << std::endl;
    mean -= controls_mean;
    controls_space.colwise() -= controls_mean;

    Rcpp::Rcerr << "SVD started" << std::endl;
    Eigen::BDCSVD<Eigen::MatrixXd> svd(controls_space, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Rcpp::Rcerr << "U matrix dims: " << svd.matrixU().rows() << " by " << svd.matrixU().cols() << std::endl;
    auto U = svd.matrixU().block(0, 0, controls_space.rows(), 10);

    Rcpp::Rcerr << "U matrix dims: " << U.rows() << " by " << U.cols() << std::endl;

    Vector rs_mean = U.transpose() * mean;
    Matrix rs_directions = U.transpose() * directions;

    std::vector<int> sel_cols;
    for (int i = 0; i < rs_directions.cols(); i++) {
        if (rs_directions.col(i).norm() > EPS) {
            sel_cols.push_back(i);
        }
    }
    Matrix rs_non_singular_vectors(rs_directions.rows(), sel_cols.size());
    for (int i = 0; i < sel_cols.size(); i++) {
        rs_non_singular_vectors.col(i) = rs_directions.col(sel_cols[i]);
    }

    Matrix rs_cov = rs_non_singular_vectors * rs_non_singular_vectors.transpose();

    Matrix rs_controls_gmatrix = U.transpose() * controls_space;
    Rcpp::Rcerr << "Subsampling started" << std::endl;
    subsampling = mvn::subsample(rs_controls_gmatrix, clustering, rs_mean, rs_cov);
    Rcpp::Rcerr << "Mahalanobis distances successfully have been calculated." << std::endl;
    subsampling.run(100000 / threads, 10, 1.0 / 3.0, threads);
}

Eigen::VectorXd matching::cluster_mean(int cluster) const {
    auto objects = clustering.elements(cluster);
    if (objects.empty()) {
        throw std::invalid_argument("Cluster is empty");
    }
    Vector v =  std::accumulate(objects.begin(), objects.end(), Vector(controls_space.rows()), [&](Vector sum, int obj) {
        return sum + controls_space.col(obj);
    });
    v /= objects.size();
    return v;
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

matching_results::matching_results(std::vector<int>&& prefix, std::vector<double>&& p_values,
                                   std::vector<double>&& lmbds, std::vector<int>&& lmbd_i,
                                   std::vector<int>&& pvals_number)
        : optimal_prefix(std::move(prefix)), pvals(std::move(p_values)),
          lambdas(std::move(lmbds)), lambda_i(std::move(lmbd_i)), pvals_num(std::move(pvals_number)) {}

}