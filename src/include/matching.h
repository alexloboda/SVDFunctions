#ifndef SRC_MATCHING_H
#define SRC_MATCHING_H

#include <vector>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <functional>
#include <stdexcept>

#include "subsample.h"
#include "lm.h"


namespace matching {

class Counts {
    int counts[3];
public:
    Counts() :counts{0, 0, 0} {}
    Counts(int hom, int het, int alt) :counts{hom, het, alt} {}

    const int& operator[] (size_t i) const {
        check_bounds(i);
        return counts[i];
    }
    int& operator[] (size_t i) {
        check_bounds(i);
        return counts[i];
    }

    int sum() const {
        return counts[0] + counts[1] + counts[2];
    }

private:
    void check_bounds(size_t i) const {
        if (i > 2) {
            throw std::domain_error("Out of bound");
        }
    }
};

class lambda_range {
    double lb;
    double ub;
public:
    lambda_range(double lb, double ub);
    bool in(double lambda);
    lambda_range();
    double distance(double lambda);
};

struct matching_results {
    std::vector<int> optimal_prefix;
    std::vector<double> pvals;
    std::vector<double> lambdas;
    std::vector<double> statistics;
    std::vector<int> lambda_i;
    std::vector<int> pvals_num;
    double optimal_lambda;

    matching_results(std::vector<int>&& optimal_prefix, std::vector<double>&& p_values, std::vector<double>&& lmbds,
                     std::vector<double>&& statistics, std::vector<int>&& lmbd_i, std::vector<int>&& pvals_number,
                     double oprimal_lambda);
};

class matching {
    static constexpr double EPS = 1e-6;

    std::vector<std::vector<int>> controls_gmatrix;
    std::shared_ptr<Eigen::MatrixXd> controls_space;

    mvn::Clustering clustering;
    mvn::subsample subsampling;

    std::function<void()> interrupts_checker;
    std::function<double(double)> qchisq;

    lambda_range hard_threshold;
    lambda_range soft_threshold;
public:
    matching(std::vector<std::vector<int>>&& controls, std::shared_ptr<Eigen::MatrixXd> space, mvn::Clustering clustering);
    void process_mvn(const Eigen::MatrixXd& directions, Eigen::VectorXd mean,
                     int threads, int start, int size_ub, int step, int iterations,
                     bool use_nystrom, size_t n_features,
                     bool use_hybrid = false, size_t rff_features = 256, uint32_t seed = 42u,
                     double ci_width_threshold = 0.5, double calibration_quantile = 0.9,
                     size_t min_calibration_samples = 8, size_t max_calibration_history = 256);
    void set_qchi_sq_function(const std::function<double(double)>& f);
    matching_results match(const std::vector<Counts>& case_counts, unsigned min_controls = 1, double min_call_rate = 0.95);

    void set_interrupts_checker(std::function<void()> checker) {
        interrupts_checker = checker;
    }

    void set_soft_threshold(lambda_range range);
    void set_hard_threshold(lambda_range range);

    size_t sa_solutions() const {
        return subsampling.solutions();
    }

    size_t sa_solution_size(size_t k) const {
        return subsampling.solution_size(k);
    }

    size_t sa_total_swaps(size_t k) const {
        return subsampling.total_swaps(k);
    }

    size_t sa_uncertain_swaps(size_t k) const {
        return subsampling.uncertain_swaps(k);
    }

    size_t sa_ci_resolved_swaps(size_t k) const {
        return subsampling.ci_resolved_swaps(k);
    }

    size_t sa_primary_resolved_swaps(size_t k) const {
        return subsampling.primary_resolved_swaps(k);
    }

    size_t sa_aux_ladder_levels() const {
        return subsampling.aux_ladder_levels();
    }

    size_t sa_aux_level_dim(size_t level) const {
        return subsampling.aux_level_dim(level);
    }

    size_t sa_aux_level_resolved_swaps(size_t k, size_t level) const {
        return subsampling.aux_level_resolved_swaps(k, level);
    }

    size_t sa_exact_unavailable_swaps(size_t k) const {
        return subsampling.exact_unavailable_swaps(k);
    }

    size_t sa_exact_evals(size_t k) const {
        return subsampling.exact_evals(k);
    }

    size_t sa_exact_evals_on_improving(size_t k) const {
        return subsampling.exact_evals_on_improving(k);
    }

    size_t sa_exact_evals_on_worsening(size_t k) const {
        return subsampling.exact_evals_on_worsening(k);
    }

    size_t sa_exact_failures(size_t k) const {
        return subsampling.exact_failures(k);
    }

    size_t sa_primary_calibration_points(size_t k) const {
        return subsampling.primary_calibration_points(k);
    }

    size_t sa_aux_calibration_points(size_t k) const {
        return subsampling.aux_calibration_points(k);
    }

    double sa_mean_primary_ci_width(size_t k) const {
        return subsampling.mean_primary_ci_width(k);
    }

    double sa_mean_aux_ci_width(size_t k) const {
        return subsampling.mean_aux_ci_width(k);
    }

    double sa_mean_selected_ci_width(size_t k) const {
        return subsampling.mean_selected_ci_width(k);
    }

    size_t sa_temperature_bins() const {
        return subsampling.temperature_bins();
    }

    size_t sa_temperature_bin_total_swaps(size_t k, size_t bin) const {
        return subsampling.temperature_bin_total_swaps(k, bin);
    }

    size_t sa_temperature_bin_accepted_swaps(size_t k, size_t bin) const {
        return subsampling.temperature_bin_accepted_swaps(k, bin);
    }

    size_t sa_temperature_bin_primary_resolved_swaps(size_t k, size_t bin) const {
        return subsampling.temperature_bin_primary_resolved_swaps(k, bin);
    }

    size_t sa_temperature_bin_aux_resolved_swaps(size_t k, size_t bin) const {
        return subsampling.temperature_bin_aux_resolved_swaps(k, bin);
    }

    size_t sa_temperature_bin_exact_evals(size_t k, size_t bin) const {
        return subsampling.temperature_bin_exact_evals(k, bin);
    }

    size_t sa_temperature_bin_exact_unavailable_swaps(size_t k, size_t bin) const {
        return subsampling.temperature_bin_exact_unavailable_swaps(k, bin);
    }

    size_t sa_temperature_bin_exact_failures(size_t k, size_t bin) const {
        return subsampling.temperature_bin_exact_failures(k, bin);
    }
private:
    double get_lambda(std::vector<double>& pvals);

    Counts count_controls(const std::vector<int>& vector, size_t j);
};

}

#endif //SRC_MATCHING_H
