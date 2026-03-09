#include "include/subsample.h"

#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include "include/third-party/cxxpool.h"

namespace mvn {

subsample::subsample(std::shared_ptr<const mvn::Matrix> X, const Clustering& clst, const mvn::Vector& mean,
                             const mvn::Matrix& cov, bool use_nystrom, size_t n_features,
                             bool use_hybrid, size_t rff_features, uint32_t seed,
                             double ci_width_threshold, double calibration_quantile,
                             size_t min_calibration_samples, size_t max_calibration_history)
            : test{std::make_shared<mvn_test>(mvn_test(X, clst, cov, mean, use_nystrom, n_features,
                                        use_hybrid, rff_features, seed))},
          ci_width_threshold(ci_width_threshold),
          calibration_quantile(calibration_quantile),
          min_calibration_samples(min_calibration_samples),
          max_calibration_history(max_calibration_history),
          clst(clst),
          wheel(std::random_device()()) {
}

namespace {

constexpr size_t TEMPERATURE_BIN_COUNT = 10;

void check_solution_vectors(const std::vector<std::vector<size_t>>& best,
                           const std::vector<double>& best_stat,
                           const std::vector<size_t>& best_size,
                           const std::vector<size_t>& total_swaps_used,
                           const std::vector<size_t>& uncertain_swaps_used,
                           const std::vector<size_t>& ci_resolved_swaps_used,
                           const std::vector<size_t>& primary_resolved_swaps_used,
                           const std::vector<std::vector<size_t>>& aux_level_resolved_swaps_used,
                           const std::vector<size_t>& exact_unavailable_swaps_used,
                           const std::vector<size_t>& exact_evals_used,
                           const std::vector<size_t>& exact_evals_on_improving_used,
                           const std::vector<size_t>& exact_evals_on_worsening_used,
                           const std::vector<size_t>& exact_eval_failures,
                           const std::vector<size_t>& primary_calibration_points_used,
                           const std::vector<size_t>& aux_calibration_points_used,
                           const std::vector<double>& mean_primary_ci_width_used,
                           const std::vector<double>& mean_aux_ci_width_used,
                           const std::vector<double>& mean_selected_ci_width_used,
                           const std::vector<std::vector<size_t>>& temperature_bin_total_swaps_used,
                           const std::vector<std::vector<size_t>>& temperature_bin_accepted_swaps_used,
                           const std::vector<std::vector<size_t>>& temperature_bin_primary_resolved_swaps_used,
                           const std::vector<std::vector<size_t>>& temperature_bin_aux_resolved_swaps_used,
                           const std::vector<std::vector<size_t>>& temperature_bin_exact_evals_used,
                           const std::vector<std::vector<size_t>>& temperature_bin_exact_unavailable_swaps_used,
                           const std::vector<std::vector<size_t>>& temperature_bin_exact_failures_used)
{
    const size_t n = best.size();
    auto require_size = [n](size_t size, const char* name) {
        if (size != n) {
            throw std::logic_error(std::string("subsample result vector '") + name +
                                   "' is inconsistent with solutions(): expected " +
                                   std::to_string(n) + ", got " + std::to_string(size));
        }
    };

    require_size(best_stat.size(), "best_stat");
    require_size(best_size.size(), "best_size");
    require_size(total_swaps_used.size(), "total_swaps_used");
    require_size(uncertain_swaps_used.size(), "uncertain_swaps_used");
    require_size(ci_resolved_swaps_used.size(), "ci_resolved_swaps_used");
    require_size(primary_resolved_swaps_used.size(), "primary_resolved_swaps_used");
    require_size(aux_level_resolved_swaps_used.size(), "aux_level_resolved_swaps_used");
    require_size(exact_unavailable_swaps_used.size(), "exact_unavailable_swaps_used");
    require_size(exact_evals_used.size(), "exact_evals_used");
    require_size(exact_evals_on_improving_used.size(), "exact_evals_on_improving_used");
    require_size(exact_evals_on_worsening_used.size(), "exact_evals_on_worsening_used");
    require_size(exact_eval_failures.size(), "exact_eval_failures");
    require_size(primary_calibration_points_used.size(), "primary_calibration_points_used");
    require_size(aux_calibration_points_used.size(), "aux_calibration_points_used");
    require_size(mean_primary_ci_width_used.size(), "mean_primary_ci_width_used");
    require_size(mean_aux_ci_width_used.size(), "mean_aux_ci_width_used");
    require_size(mean_selected_ci_width_used.size(), "mean_selected_ci_width_used");
    require_size(temperature_bin_total_swaps_used.size(), "temperature_bin_total_swaps_used");
    require_size(temperature_bin_accepted_swaps_used.size(), "temperature_bin_accepted_swaps_used");
    require_size(temperature_bin_primary_resolved_swaps_used.size(), "temperature_bin_primary_resolved_swaps_used");
    require_size(temperature_bin_aux_resolved_swaps_used.size(), "temperature_bin_aux_resolved_swaps_used");
    require_size(temperature_bin_exact_evals_used.size(), "temperature_bin_exact_evals_used");
    require_size(temperature_bin_exact_unavailable_swaps_used.size(), "temperature_bin_exact_unavailable_swaps_used");
    require_size(temperature_bin_exact_failures_used.size(), "temperature_bin_exact_failures_used");

    for (size_t i = 0; i < n; ++i) {
        for (const auto* bins : {&temperature_bin_total_swaps_used, &temperature_bin_accepted_swaps_used,
                                 &temperature_bin_primary_resolved_swaps_used, &temperature_bin_aux_resolved_swaps_used,
                                 &temperature_bin_exact_evals_used, &temperature_bin_exact_unavailable_swaps_used,
                                 &temperature_bin_exact_failures_used}) {
            if (bins->at(i).size() != TEMPERATURE_BIN_COUNT) {
                throw std::logic_error("temperature-bin diagnostics have inconsistent width");
            }
        }
    }
}

size_t temperature_bin_index(size_t iteration, size_t iterations)
{
    if (iterations == 0) {
        return 0;
    }
    return std::min(TEMPERATURE_BIN_COUNT - 1, (iteration * TEMPERATURE_BIN_COUNT) / iterations);
}

struct run_result {
    std::shared_ptr<mvn_test> test;
    size_t total_swaps = 0;
    size_t uncertain_swaps = 0;
    size_t ci_resolved_swaps = 0;
    size_t primary_resolved_swaps = 0;
    std::vector<size_t> aux_level_resolved_swaps;
    size_t exact_unavailable_swaps = 0;
    size_t exact_evals = 0;
    size_t exact_evals_on_improving = 0;
    size_t exact_evals_on_worsening = 0;
    size_t exact_failures = 0;
    size_t primary_calibration_points = 0;
    size_t aux_calibration_points = 0;
    double primary_ci_width_sum = 0.0;
    double aux_ci_width_sum = 0.0;
    double selected_ci_width_sum = 0.0;
    size_t ci_width_observations = 0;
    std::vector<size_t> temperature_bin_total_swaps;
    std::vector<size_t> temperature_bin_accepted_swaps;
    std::vector<size_t> temperature_bin_primary_resolved_swaps;
    std::vector<size_t> temperature_bin_aux_resolved_swaps;
    std::vector<size_t> temperature_bin_exact_evals;
    std::vector<size_t> temperature_bin_exact_unavailable_swaps;
    std::vector<size_t> temperature_bin_exact_failures;
};

template <class T>
constexpr std::add_const_t<T>& const_ref(T& t) noexcept
{
    return t;
}

double acceptance_probability(double delta, double temperature)
{
    if (delta <= 0.0) {
        return 1.0;
    }
    return std::exp(-delta / temperature);
}

double empirical_quantile(const std::vector<double>& values, double probability)
{
    if (values.empty()) {
        return 1.0;
    }

    std::vector<double> sorted(values.begin(), values.end());
    std::sort(sorted.begin(), sorted.end());

    const double clamped_probability = std::min(1.0, std::max(0.0, probability));
    const double position = clamped_probability * static_cast<double>(sorted.size() - 1);
    const size_t lower = static_cast<size_t>(std::floor(position));
    const size_t upper = static_cast<size_t>(std::ceil(position));
    if (lower == upper) {
        return sorted[lower];
    }

    const double weight = position - static_cast<double>(lower);
    return sorted[lower] + weight * (sorted[upper] - sorted[lower]);
}

double calibration_half_width(const std::vector<double>& residuals, double calibration_quantile,
                              size_t min_calibration_samples)
{
    if (residuals.size() < min_calibration_samples) {
        return 1.0;
    }
    return empirical_quantile(residuals, calibration_quantile);
}

double interval_width(double center, double half_width)
{
    const double lower = std::max(0.0, center - half_width);
    const double upper = std::min(1.0, center + half_width);
    return upper - lower;
}

void push_bounded(std::vector<double>& values, double value, size_t max_calibration_history)
{
    if (values.size() >= max_calibration_history) {
        values.erase(values.begin());
    }
    values.push_back(value);
}

}

void subsample::run(size_t iterations, size_t restarts, double t_start, double c, size_t pool_size, size_t start, size_t size_ub,
                    size_t step) {
    using std::vector;
    using namespace std::chrono_literals;

    best.clear();
    best_stat.clear();
    best_size.clear();
    total_swaps_used.clear();
    uncertain_swaps_used.clear();
    ci_resolved_swaps_used.clear();
    primary_resolved_swaps_used.clear();
    aux_level_resolved_swaps_used.clear();
    exact_unavailable_swaps_used.clear();
    exact_evals_used.clear();
    exact_evals_on_improving_used.clear();
    exact_evals_on_worsening_used.clear();
    exact_eval_failures.clear();
    primary_calibration_points_used.clear();
    aux_calibration_points_used.clear();
    mean_primary_ci_width_used.clear();
    mean_aux_ci_width_used.clear();
    mean_selected_ci_width_used.clear();
    temperature_bin_total_swaps_used.clear();
    temperature_bin_accepted_swaps_used.clear();
    temperature_bin_primary_resolved_swaps_used.clear();
    temperature_bin_aux_resolved_swaps_used.clear();
    temperature_bin_exact_evals_used.clear();
    temperature_bin_exact_unavailable_swaps_used.clear();
    temperature_bin_exact_failures_used.clear();

    cxxpool::thread_pool pool(pool_size);
    size_t curr_size = start;

    while (curr_size <= size_ub) {
        Rcpp::checkUserInterrupt();
        std::vector<std::future<run_result>> thread_solutions;
        for (size_t t = 0; t < restarts; t++) {
            thread_solutions.push_back(pool.push([test = const_ref(test), iterations, t_start, curr_size, c,
                                                         seed = wheel(),
                                                         ci_width_threshold = this->ci_width_threshold,
                                                         calibration_quantile = this->calibration_quantile,
                                                         min_calibration_samples = this->min_calibration_samples,
                                                         max_calibration_history = this->max_calibration_history]() -> run_result {
                run_result result;
                double t = t_start;
                std::mt19937 mersenne_wheel(seed);
                std::shared_ptr<mvn_test> local_test = test->clone();
                std::vector<double> primary_residuals;
                std::vector<std::vector<double>> aux_residuals(local_test->aux_statistic_levels());
                result.aux_level_resolved_swaps.assign(local_test->aux_statistic_levels(), 0);
                result.temperature_bin_total_swaps.assign(TEMPERATURE_BIN_COUNT, 0);
                result.temperature_bin_accepted_swaps.assign(TEMPERATURE_BIN_COUNT, 0);
                result.temperature_bin_primary_resolved_swaps.assign(TEMPERATURE_BIN_COUNT, 0);
                result.temperature_bin_aux_resolved_swaps.assign(TEMPERATURE_BIN_COUNT, 0);
                result.temperature_bin_exact_evals.assign(TEMPERATURE_BIN_COUNT, 0);
                result.temperature_bin_exact_unavailable_swaps.assign(TEMPERATURE_BIN_COUNT, 0);
                result.temperature_bin_exact_failures.assign(TEMPERATURE_BIN_COUNT, 0);
                result.total_swaps = iterations;
                while (local_test->subsample_size() < curr_size) {
                    local_test->add_one();
                }
                std::uniform_real_distribution<double> random_unif(0.0, 1.0);
                double score = local_test->get_normality_statistic();
                std::vector<double> aux_scores(local_test->aux_statistic_levels(), score);
                if (local_test->has_aux_statistic()) {
                    for (size_t level = 0; level < aux_scores.size(); ++level) {
                        aux_scores[level] = local_test->get_aux_normality_statistic(level);
                    }
                }
                for (size_t k = 0; k < iterations; k++) {
                    t = c * t;
                    const size_t temp_bin = temperature_bin_index(k, iterations);
                    result.temperature_bin_total_swaps[temp_bin]++;
                    local_test->swap_once();
                    double new_score = local_test->get_normality_statistic();
                    std::vector<double> new_aux_scores(aux_scores.size(), new_score);
                    if (local_test->has_aux_statistic()) {
                        for (size_t level = 0; level < new_aux_scores.size(); ++level) {
                            new_aux_scores[level] = local_test->get_aux_normality_statistic(level);
                        }
                    }

                    const double delta_primary = new_score - score;
                    const bool primary_improving = delta_primary <= 0.0;
                    const double p_primary = acceptance_probability(delta_primary, t);

                    bool accept = true;

                    if (local_test->has_aux_statistic()) {
                        std::vector<double> p_aux_values(new_aux_scores.size(), p_primary);
                        for (size_t level = 0; level < new_aux_scores.size(); ++level) {
                            const double delta_aux = new_aux_scores[level] - aux_scores[level];
                            p_aux_values[level] = acceptance_probability(delta_aux, t);
                        }

                        const double primary_local_spread = std::abs(p_primary - p_aux_values.front());
                        const double primary_half_width = std::max(
                            primary_local_spread,
                            calibration_half_width(primary_residuals, calibration_quantile, min_calibration_samples));
                        const double primary_width = interval_width(p_primary, primary_half_width);

                        double representative_aux_width = 1.0;
                        double selected_width = primary_width;
                        double cheap_probability = p_primary;
                        double fallback_width = primary_width;
                        double fallback_probability = p_primary;
                        bool resolved_by_ci = (primary_width <= ci_width_threshold);
                        bool resolved_by_primary = resolved_by_ci;
                        bool resolved_by_aux = false;
                        size_t resolved_aux_level = p_aux_values.size();

                        double reference_probability = p_primary;
                        for (size_t level = 0; level < p_aux_values.size(); ++level) {
                            const double aux_half_width = std::max(
                                std::abs(p_aux_values[level] - reference_probability),
                                calibration_half_width(aux_residuals[level], calibration_quantile, min_calibration_samples));
                            const double aux_width = interval_width(p_aux_values[level], aux_half_width);
                            representative_aux_width = aux_width;
                            if (aux_width < fallback_width) {
                                fallback_width = aux_width;
                                fallback_probability = p_aux_values[level];
                            }
                            if (aux_width <= ci_width_threshold) {
                                if (!resolved_by_ci || aux_width < selected_width) {
                                    selected_width = aux_width;
                                    cheap_probability = p_aux_values[level];
                                }
                                resolved_by_ci = true;
                                resolved_by_primary = false;
                                resolved_by_aux = true;
                                resolved_aux_level = level;
                                break;
                            }
                            reference_probability = p_aux_values[level];
                        }

                        result.primary_ci_width_sum += primary_width;
                        result.aux_ci_width_sum += representative_aux_width;
                        result.ci_width_observations++;

                        if (!resolved_by_ci) {
                            selected_width = fallback_width;
                            cheap_probability = fallback_probability;
                        }
                        result.selected_ci_width_sum += selected_width;

                        const bool need_exact = !resolved_by_ci;

                        if (need_exact) {
                            result.uncertain_swaps++;
                            if (local_test->last_swap_has_equal_effect_size()) {
                                try {
                                    result.exact_evals++;
                                    result.temperature_bin_exact_evals[temp_bin]++;
                                    if (primary_improving) {
                                        result.exact_evals_on_improving++;
                                    } else {
                                        result.exact_evals_on_worsening++;
                                    }

                                    const double delta_exact = local_test->exact_delta_last_swap();
                                    const double p_exact = acceptance_probability(delta_exact, t);
                                    push_bounded(primary_residuals, std::abs(p_exact - p_primary), max_calibration_history);
                                    for (size_t level = 0; level < p_aux_values.size(); ++level) {
                                        push_bounded(aux_residuals[level], std::abs(p_exact - p_aux_values[level]), max_calibration_history);
                                    }
                                    result.primary_calibration_points = primary_residuals.size();
                                    result.aux_calibration_points = aux_residuals.empty() ? 0 : aux_residuals.back().size();
                                    accept = (random_unif(mersenne_wheel) < p_exact);
                                } catch (...) {
                                    result.exact_failures++;
                                    result.temperature_bin_exact_failures[temp_bin]++;
                                    accept = (random_unif(mersenne_wheel) < cheap_probability);
                                }
                            } else {
                                result.exact_unavailable_swaps++;
                                result.temperature_bin_exact_unavailable_swaps[temp_bin]++;
                                accept = (random_unif(mersenne_wheel) < cheap_probability);
                            }
                        } else {
                            result.ci_resolved_swaps++;
                            if (resolved_by_primary) {
                                result.primary_resolved_swaps++;
                                result.temperature_bin_primary_resolved_swaps[temp_bin]++;
                            }
                            if (resolved_by_aux && resolved_aux_level < result.aux_level_resolved_swaps.size()) {
                                result.aux_level_resolved_swaps[resolved_aux_level]++;
                                result.temperature_bin_aux_resolved_swaps[temp_bin]++;
                            }
                            accept = (random_unif(mersenne_wheel) < cheap_probability);
                        }
                    } else {
                        accept = (random_unif(mersenne_wheel) < p_primary);
                    }

                    if (!accept) {
                        local_test->swap_once(true);
                    } else {
                        result.temperature_bin_accepted_swaps[temp_bin]++;
                        score = new_score;
                        aux_scores = new_aux_scores;
                    }
                }
                result.test = local_test;
                return result;
            }));
        }

        size_t total_swaps_curr = 0;
        size_t uncertain_swaps_curr = 0;
        size_t ci_resolved_swaps_curr = 0;
        size_t primary_resolved_swaps_curr = 0;
        std::vector<size_t> aux_level_resolved_swaps_curr(test->aux_statistic_levels(), 0);
        size_t exact_unavailable_swaps_curr = 0;
        size_t exact_evals_curr = 0;
        size_t exact_evals_improving_curr = 0;
        size_t exact_evals_worsening_curr = 0;
        size_t exact_failures_curr = 0;
        size_t primary_calibration_points_curr = 0;
        size_t aux_calibration_points_curr = 0;
        double primary_ci_width_sum_curr = 0.0;
        double aux_ci_width_sum_curr = 0.0;
        double selected_ci_width_sum_curr = 0.0;
        size_t ci_width_observations_curr = 0;
        std::vector<size_t> temperature_bin_total_swaps_curr(TEMPERATURE_BIN_COUNT, 0);
        std::vector<size_t> temperature_bin_accepted_swaps_curr(TEMPERATURE_BIN_COUNT, 0);
        std::vector<size_t> temperature_bin_primary_resolved_swaps_curr(TEMPERATURE_BIN_COUNT, 0);
        std::vector<size_t> temperature_bin_aux_resolved_swaps_curr(TEMPERATURE_BIN_COUNT, 0);
        std::vector<size_t> temperature_bin_exact_evals_curr(TEMPERATURE_BIN_COUNT, 0);
        std::vector<size_t> temperature_bin_exact_unavailable_swaps_curr(TEMPERATURE_BIN_COUNT, 0);
        std::vector<size_t> temperature_bin_exact_failures_curr(TEMPERATURE_BIN_COUNT, 0);

        for (size_t i = 0; i < thread_solutions.size(); i++) {
            auto& future = thread_solutions[i];
            run_result run = future.get();
            std::shared_ptr<mvn_test> thread = run.test;
            total_swaps_curr += run.total_swaps;
            uncertain_swaps_curr += run.uncertain_swaps;
            ci_resolved_swaps_curr += run.ci_resolved_swaps;
            primary_resolved_swaps_curr += run.primary_resolved_swaps;
            if (aux_level_resolved_swaps_curr.size() < run.aux_level_resolved_swaps.size()) {
                aux_level_resolved_swaps_curr.resize(run.aux_level_resolved_swaps.size(), 0);
            }
            for (size_t level = 0; level < run.aux_level_resolved_swaps.size(); ++level) {
                aux_level_resolved_swaps_curr[level] += run.aux_level_resolved_swaps[level];
            }
            exact_unavailable_swaps_curr += run.exact_unavailable_swaps;
            exact_evals_curr += run.exact_evals;
            exact_evals_improving_curr += run.exact_evals_on_improving;
            exact_evals_worsening_curr += run.exact_evals_on_worsening;
            exact_failures_curr += run.exact_failures;
            primary_calibration_points_curr += run.primary_calibration_points;
            aux_calibration_points_curr += run.aux_calibration_points;
            primary_ci_width_sum_curr += run.primary_ci_width_sum;
            aux_ci_width_sum_curr += run.aux_ci_width_sum;
            selected_ci_width_sum_curr += run.selected_ci_width_sum;
            ci_width_observations_curr += run.ci_width_observations;
            for (size_t bin = 0; bin < TEMPERATURE_BIN_COUNT; ++bin) {
                temperature_bin_total_swaps_curr[bin] += run.temperature_bin_total_swaps[bin];
                temperature_bin_accepted_swaps_curr[bin] += run.temperature_bin_accepted_swaps[bin];
                temperature_bin_primary_resolved_swaps_curr[bin] += run.temperature_bin_primary_resolved_swaps[bin];
                temperature_bin_aux_resolved_swaps_curr[bin] += run.temperature_bin_aux_resolved_swaps[bin];
                temperature_bin_exact_evals_curr[bin] += run.temperature_bin_exact_evals[bin];
                temperature_bin_exact_unavailable_swaps_curr[bin] += run.temperature_bin_exact_unavailable_swaps[bin];
                temperature_bin_exact_failures_curr[bin] += run.temperature_bin_exact_failures[bin];
            }
            if (i == 0) {
                test = thread;
            } else if (*thread < *test) {
                test = thread;
            }
        }

        best.push_back(test->current_subset());
        best_stat.push_back(test->get_normality_statistic());
        best_size.push_back(curr_size);
        total_swaps_used.push_back(total_swaps_curr);
        uncertain_swaps_used.push_back(uncertain_swaps_curr);
        ci_resolved_swaps_used.push_back(ci_resolved_swaps_curr);
        primary_resolved_swaps_used.push_back(primary_resolved_swaps_curr);
        aux_level_resolved_swaps_used.push_back(aux_level_resolved_swaps_curr);
        exact_unavailable_swaps_used.push_back(exact_unavailable_swaps_curr);
        exact_evals_used.push_back(exact_evals_curr);
        exact_evals_on_improving_used.push_back(exact_evals_improving_curr);
        exact_evals_on_worsening_used.push_back(exact_evals_worsening_curr);
        exact_eval_failures.push_back(exact_failures_curr);
        primary_calibration_points_used.push_back(primary_calibration_points_curr);
        aux_calibration_points_used.push_back(aux_calibration_points_curr);
        if (ci_width_observations_curr == 0) {
            mean_primary_ci_width_used.push_back(NA_REAL);
            mean_aux_ci_width_used.push_back(NA_REAL);
            mean_selected_ci_width_used.push_back(NA_REAL);
        } else {
            const double denominator = static_cast<double>(ci_width_observations_curr);
            mean_primary_ci_width_used.push_back(primary_ci_width_sum_curr / denominator);
            mean_aux_ci_width_used.push_back(aux_ci_width_sum_curr / denominator);
            mean_selected_ci_width_used.push_back(selected_ci_width_sum_curr / denominator);
        }
        temperature_bin_total_swaps_used.push_back(temperature_bin_total_swaps_curr);
        temperature_bin_accepted_swaps_used.push_back(temperature_bin_accepted_swaps_curr);
        temperature_bin_primary_resolved_swaps_used.push_back(temperature_bin_primary_resolved_swaps_curr);
        temperature_bin_aux_resolved_swaps_used.push_back(temperature_bin_aux_resolved_swaps_curr);
        temperature_bin_exact_evals_used.push_back(temperature_bin_exact_evals_curr);
        temperature_bin_exact_unavailable_swaps_used.push_back(temperature_bin_exact_unavailable_swaps_curr);
        temperature_bin_exact_failures_used.push_back(temperature_bin_exact_failures_curr);

        curr_size += step;

        if (curr_size > test->sample_size()) {
            return;
        }
    }
}

size_t subsample::solutions() const {
    check_solution_vectors(best, best_stat, best_size, total_swaps_used, uncertain_swaps_used,
                           ci_resolved_swaps_used, primary_resolved_swaps_used, aux_level_resolved_swaps_used,
                           exact_unavailable_swaps_used, exact_evals_used,
                           exact_evals_on_improving_used, exact_evals_on_worsening_used,
                           exact_eval_failures, primary_calibration_points_used,
                           aux_calibration_points_used, mean_primary_ci_width_used,
                           mean_aux_ci_width_used, mean_selected_ci_width_used,
                           temperature_bin_total_swaps_used, temperature_bin_accepted_swaps_used,
                           temperature_bin_primary_resolved_swaps_used, temperature_bin_aux_resolved_swaps_used,
                           temperature_bin_exact_evals_used, temperature_bin_exact_unavailable_swaps_used,
                           temperature_bin_exact_failures_used);
    return best.size();
}

std::vector<size_t> subsample::get_solution(size_t k) const {
    return best.at(k);
}

subsample::subsample() {}

double subsample::statistic(size_t k) {
    return best_stat.at(k);
}

size_t subsample::solution_size(size_t k) const {
    return best_size.at(k);
}

size_t subsample::total_swaps(size_t k) const {
    return total_swaps_used.at(k);
}

size_t subsample::uncertain_swaps(size_t k) const {
    return uncertain_swaps_used.at(k);
}

size_t subsample::ci_resolved_swaps(size_t k) const {
    return ci_resolved_swaps_used.at(k);
}

size_t subsample::primary_resolved_swaps(size_t k) const {
    return primary_resolved_swaps_used.at(k);
}

size_t subsample::aux_ladder_levels() const {
    return aux_level_resolved_swaps_used.empty() ? 0 : aux_level_resolved_swaps_used.front().size();
}

size_t subsample::aux_level_resolved_swaps(size_t k, size_t level) const {
    return aux_level_resolved_swaps_used.at(k).at(level);
}

size_t subsample::exact_unavailable_swaps(size_t k) const {
    return exact_unavailable_swaps_used.at(k);
}

size_t subsample::exact_evals(size_t k) const {
    return exact_evals_used.at(k);
}

size_t subsample::exact_evals_on_improving(size_t k) const {
    return exact_evals_on_improving_used.at(k);
}

size_t subsample::exact_evals_on_worsening(size_t k) const {
    return exact_evals_on_worsening_used.at(k);
}

size_t subsample::exact_failures(size_t k) const {
    return exact_eval_failures.at(k);
}

size_t subsample::primary_calibration_points(size_t k) const {
    return primary_calibration_points_used.at(k);
}

size_t subsample::aux_calibration_points(size_t k) const {
    return aux_calibration_points_used.at(k);
}

double subsample::mean_primary_ci_width(size_t k) const {
    return mean_primary_ci_width_used.at(k);
}

double subsample::mean_aux_ci_width(size_t k) const {
    return mean_aux_ci_width_used.at(k);
}

double subsample::mean_selected_ci_width(size_t k) const {
    return mean_selected_ci_width_used.at(k);
}

size_t subsample::temperature_bins() const {
    return TEMPERATURE_BIN_COUNT;
}

size_t subsample::temperature_bin_total_swaps(size_t k, size_t bin) const {
    return temperature_bin_total_swaps_used.at(k).at(bin);
}

size_t subsample::temperature_bin_accepted_swaps(size_t k, size_t bin) const {
    return temperature_bin_accepted_swaps_used.at(k).at(bin);
}

size_t subsample::temperature_bin_primary_resolved_swaps(size_t k, size_t bin) const {
    return temperature_bin_primary_resolved_swaps_used.at(k).at(bin);
}

size_t subsample::temperature_bin_aux_resolved_swaps(size_t k, size_t bin) const {
    return temperature_bin_aux_resolved_swaps_used.at(k).at(bin);
}

size_t subsample::temperature_bin_exact_evals(size_t k, size_t bin) const {
    return temperature_bin_exact_evals_used.at(k).at(bin);
}

size_t subsample::temperature_bin_exact_unavailable_swaps(size_t k, size_t bin) const {
    return temperature_bin_exact_unavailable_swaps_used.at(k).at(bin);
}

size_t subsample::temperature_bin_exact_failures(size_t k, size_t bin) const {
    return temperature_bin_exact_failures_used.at(k).at(bin);
}

}