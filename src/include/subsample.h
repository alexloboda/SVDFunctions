#ifndef SRC_SUBSAMPLE_H
#define SRC_SUBSAMPLE_H

#include <vector>
#include "mvn_test.h"
namespace mvn {

class subsample {
    std::shared_ptr<mvn_test> test;

    double ci_width_threshold = 0.5;
    double calibration_quantile = 0.9;
    size_t min_calibration_samples = 8;
    size_t max_calibration_history = 256;

    std::vector<std::vector<size_t>> best;
    std::vector<double> best_stat;
    std::vector<size_t> best_size;
    std::vector<size_t> total_swaps_used;
    std::vector<size_t> uncertain_swaps_used;
    std::vector<size_t> ci_resolved_swaps_used;
    std::vector<size_t> exact_unavailable_swaps_used;
    std::vector<size_t> exact_evals_used;
    std::vector<size_t> exact_evals_on_improving_used;
    std::vector<size_t> exact_evals_on_worsening_used;
    std::vector<size_t> exact_eval_failures;
    std::vector<size_t> primary_calibration_points_used;
    std::vector<size_t> aux_calibration_points_used;
    std::vector<double> mean_primary_ci_width_used;
    std::vector<double> mean_aux_ci_width_used;
    std::vector<double> mean_selected_ci_width_used;
    Clustering clst;

    std::mt19937 wheel;
public:
    subsample();
    subsample(std::shared_ptr<const Matrix> X, const Clustering& clst, const Vector& mean, const Matrix& cov,
              bool use_nystrom, size_t n_features,
              bool use_hybrid = false, size_t rff_features = 256, uint32_t seed = 42u,
              double ci_width_threshold = 0.5, double calibration_quantile = 0.9,
              size_t min_calibration_samples = 8, size_t max_calibration_history = 256);
    subsample(subsample&&) = default;
    subsample& operator=(subsample&& other) = default;

    void run(size_t iterations, size_t restarts, double t0, double c, size_t threads, size_t start, size_t size_ub,
             size_t step);

    std::vector<size_t> get_solution(size_t size) const;
    size_t solutions() const;
    double statistic(size_t k);
    size_t solution_size(size_t k) const;
    size_t total_swaps(size_t k) const;
    size_t uncertain_swaps(size_t k) const;
    size_t ci_resolved_swaps(size_t k) const;
    size_t exact_unavailable_swaps(size_t k) const;
    size_t exact_evals(size_t k) const;
    size_t exact_evals_on_improving(size_t k) const;
    size_t exact_evals_on_worsening(size_t k) const;
    size_t exact_failures(size_t k) const;
    size_t primary_calibration_points(size_t k) const;
    size_t aux_calibration_points(size_t k) const;
    double mean_primary_ci_width(size_t k) const;
    double mean_aux_ci_width(size_t k) const;
    double mean_selected_ci_width(size_t k) const;
};

}
#endif //SRC_SUBSAMPLE_H
