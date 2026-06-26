#include "include/matching.h"
#include "include/hw.h"
#include "include/lm.h"

#include <numeric>
#include <utility>
#include <fstream>
#include <Rcpp.h>
#include <cmath>

namespace Eigen{
template<class Matrix>
void write_binary(const char* filename, const Matrix& matrix){
    std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
    typename Matrix::Index rows=matrix.rows(), cols=matrix.cols();
    out.write((char*) (&rows), sizeof(typename Matrix::Index));
    out.write((char*) (&cols), sizeof(typename Matrix::Index));
    out.write((char*) matrix.data(), rows*cols*sizeof(typename Matrix::Scalar) );
    out.close();
}
template<class Matrix>
void read_binary(const char* filename, Matrix& matrix){
    std::ifstream in(filename, std::ios::in | std::ios::binary);
    typename Matrix::Index rows=0, cols=0;
    in.read((char*) (&rows),sizeof(typename Matrix::Index));
    in.read((char*) (&cols),sizeof(typename Matrix::Index));
    matrix.resize(rows, cols);
    in.read( (char *) matrix.data() , rows*cols*sizeof(typename Matrix::Scalar) );
    in.close();
}
} // Eigen::

namespace {

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

}

namespace matching {

matching::matching(std::vector<std::vector<ClusterCounts>>&& cluster_counts,
                   std::shared_ptr<Eigen::MatrixXd> controls_space,
                   mvn::Clustering clustering) : cluster_counts(std::move(cluster_counts)),
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

matching_results matching::match(const std::vector<Counts>& case_counts, unsigned min_controls, double min_call_rate) {
    size_t n_variants = case_counts.size();
    auto control_sample_size = [this](const std::vector<size_t>& groups) {
        size_t total = 0;
        for (size_t group: groups) {
            total += clustering.cluster_size(group);
        }
        return total;
    };

    std::vector<bool> snp_mask = check_user_counts(case_counts);
    auto lms = init_lms(case_counts);

    double lambda = std::numeric_limits<double>::infinity();

    std::vector<double> lambdas;
    std::vector<double> stats;

    std::vector<int> optimal_clusters;
    std::vector<size_t> optimal_groups;
    std::vector<int> lambda_i;
    std::vector<int> pvals_num;
    std::vector<double> optimal_pvals;

    for (size_t k = 0, step = 0; k < subsampling.solutions(); k++, step++) {
        if (step % 100 == 0) {
            interrupts_checker();
        }

        stats.push_back(subsampling.statistic(k));

        std::vector<double> pvals;

        auto control_groups = subsampling.get_solution(k);
        if (control_groups.size() < min_controls) {
            continue;
        }
        size_t controls_size = control_sample_size(control_groups);

        for (size_t j = 0; j < n_variants; j++) {
            if (!snp_mask[j]) {
                continue;
            }

            Counts controls_counts = count_controls(control_groups, j);

            auto overall = controls_counts[0] + controls_counts[1] + controls_counts[2];
            if ((double)overall / (double)controls_size < min_call_rate) {
                continue;
            }

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

        if (controls_size >= min_controls && pvals.size() > 10) {
            double cur_lambda = get_lambda(pvals);

            lambdas.push_back(cur_lambda);
            lambda_i.push_back(static_cast<int>(controls_size));
            pvals_num.push_back(pvals.size());

            if (hard_threshold.in(cur_lambda)) {
                if (soft_threshold.in(cur_lambda) ||
                        (!soft_threshold.in(lambda) &&
                        soft_threshold.distance(cur_lambda) < soft_threshold.distance(lambda))) {
                    optimal_groups = control_groups;
                    lambda = cur_lambda;
                    optimal_pvals = pvals;
                }
            }
        }
    }

    if (!optimal_groups.empty()) {
        optimal_clusters.reserve(optimal_groups.size());
        for (size_t group: optimal_groups) {
            optimal_clusters.push_back(static_cast<int>(group));
        }
    }

    return {std::move(optimal_clusters), std::move(optimal_pvals), std::move(lambdas),
            std::move(stats), std::move(lambda_i), std::move(pvals_num), lambda};
}

void matching::process_mvn(const Matrix& directions, Vector mean,
                           int sa_threads, int start, int ub, int step, int iterations,
                           int exact_precompute_threads, int exact_cluster_tile_size) {
    const double EPS = 1e-18;
    Rcpp::Rcerr << "Starting processing controls space." << std::endl;
    Rcpp::Rcerr << "The size of controls space is " << controls_space->rows() << " by " << controls_space->cols() << std::endl;

    mvn::PrecomputeConfig config;
    if (exact_precompute_threads > 0) {
        config.threads = static_cast<size_t>(exact_precompute_threads);
    }
    if (exact_cluster_tile_size > 0) {
        config.cluster_tile_size = static_cast<size_t>(exact_cluster_tile_size);
    }

    {
        Matrix rs_cov = directions * directions.transpose();
        subsampling = mvn::subsample(controls_space, clustering, mean, rs_cov, config);
    }
    controls_space.reset();
    Rcpp::Rcerr << "Mahalanobis distances have been successfully calculated." << std::endl;
    double c = std::pow(EPS, 1.0 / (double)iterations);
    subsampling.run(iterations, 4, 1.0 , c, sa_threads, start, ub, step);
}

Counts matching::count_controls(const std::vector<size_t>& groups, size_t variant) {
    Counts counts;
    for (size_t group: groups) {
        const ClusterCounts& cluster = cluster_counts[group][variant];
        counts[0] += cluster[0];
        counts[1] += cluster[1];
        counts[2] += cluster[2];
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
                                   std::vector<double>&& lmbds, std::vector<double>&& stats,
                                   std::vector<int>&& lmbd_i,
                                   std::vector<int>&& pvals_number, double optimal_lambda)
        : optimal_prefix(std::move(prefix)), pvals(std::move(p_values)),
          lambdas(std::move(lmbds)), statistics(std::move(stats)), lambda_i(std::move(lmbd_i)),
          pvals_num(std::move(pvals_number)), optimal_lambda(optimal_lambda) {}

}