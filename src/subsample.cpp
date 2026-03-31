#include "include/subsample.h"

#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <limits>
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
constexpr double HYBRID_RESOLVED_AUDIT_RATE = 0.05;
constexpr double HYBRID_SHADOW_AUDIT_RATE = 0.05;
constexpr double ADAPTIVE_BIN_CENTER_ABS_QUANTILE = 0.35;
constexpr double ADAPTIVE_BIN_SIDE_SPLIT_QUANTILE = 0.5;
constexpr size_t CALIBRATION_LOOKUP_SOURCE_COUNT = 4;

enum class calibration_lookup_source {
    local,
    side,
    global,
    insufficient
};

struct calibration_lookup_result {
    double half_width = std::numeric_limits<double>::infinity();
    size_t sample_count = 0;
    calibration_lookup_source source = calibration_lookup_source::insufficient;
};

struct calibration_record {
    double cheap_delta = 0.0;
    double abs_error = 0.0;
};

struct adaptive_residual_bins {
    size_t max_history = 256;
    std::vector<calibration_record> records;

    adaptive_residual_bins() = default;

    explicit adaptive_residual_bins(size_t max_history)
        : max_history(max_history) {}

    void observe(double cheap_delta, double exact_delta) {
        if (records.size() >= max_history) {
            records.erase(records.begin());
        }
        records.push_back({cheap_delta, std::abs(exact_delta - cheap_delta)});
    }

    size_t size() const {
        return records.size();
    }

    calibration_lookup_result half_width(double cheap_delta, double calibration_quantile,
                                         size_t min_calibration_samples) const;
};

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
                           const std::vector<size_t>& primary_local_calibration_lookups_used,
                           const std::vector<size_t>& primary_side_calibration_lookups_used,
                           const std::vector<size_t>& primary_global_calibration_lookups_used,
                           const std::vector<size_t>& primary_insufficient_calibration_lookups_used,
                           const std::vector<size_t>& aux_local_calibration_lookups_used,
                           const std::vector<size_t>& aux_side_calibration_lookups_used,
                           const std::vector<size_t>& aux_global_calibration_lookups_used,
                           const std::vector<size_t>& aux_insufficient_calibration_lookups_used,
                           const std::vector<double>& mean_primary_ci_width_used,
                           const std::vector<double>& mean_aux_ci_width_used,
                           const std::vector<double>& mean_selected_ci_width_used,
                           const std::vector<double>& mean_scanned_aux_levels_used,
                           const std::vector<size_t>& shadow_audits_used,
                           const std::vector<size_t>& shadow_exact_failures_used,
                           const std::vector<size_t>& shadow_primary_decision_mismatches_used,
                           const std::vector<size_t>& shadow_selected_decision_mismatches_used,
                           const std::vector<size_t>& shadow_selected_interval_hits_used,
                           const std::vector<size_t>& shadow_full_scan_better_swaps_used,
                           const std::vector<double>& mean_shadow_selected_regret_used,
                           const std::vector<double>& mean_shadow_primary_abs_delta_error_used,
                           const std::vector<double>& mean_shadow_primary_delta_bias_used,
                           const std::vector<double>& mean_shadow_primary_p_bias_used,
                           const std::vector<double>& mean_shadow_selected_abs_delta_error_used,
                           const std::vector<double>& mean_shadow_selected_abs_p_error_used,
                           const std::vector<double>& mean_shadow_selected_delta_bias_used,
                           const std::vector<double>& mean_shadow_selected_p_bias_used,
                           const std::vector<std::vector<size_t>>& shadow_selected_source_audits_used,
                           const std::vector<std::vector<size_t>>& shadow_selected_source_exact_failures_used,
                           const std::vector<std::vector<size_t>>& shadow_selected_source_decision_mismatches_used,
                           const std::vector<std::vector<size_t>>& shadow_selected_source_interval_hits_used,
                           const std::vector<std::vector<size_t>>& shadow_selected_source_full_scan_better_swaps_used,
                           const std::vector<std::vector<double>>& mean_shadow_selected_source_regret_used,
                           const std::vector<std::vector<double>>& mean_shadow_selected_source_abs_delta_error_used,
                           const std::vector<std::vector<double>>& mean_shadow_selected_source_abs_p_error_used,
                           const std::vector<std::vector<double>>& mean_shadow_selected_source_delta_bias_used,
                           const std::vector<std::vector<double>>& mean_shadow_selected_source_p_bias_used,
                           const std::vector<std::vector<size_t>>& temperature_bin_total_swaps_used,
                           const std::vector<std::vector<size_t>>& temperature_bin_accepted_swaps_used,
                           const std::vector<std::vector<size_t>>& temperature_bin_primary_resolved_swaps_used,
                           const std::vector<std::vector<size_t>>& temperature_bin_aux_resolved_swaps_used,
                           const std::vector<std::vector<size_t>>& temperature_bin_exact_evals_used,
                           const std::vector<std::vector<size_t>>& temperature_bin_exact_unavailable_swaps_used,
                           const std::vector<std::vector<size_t>>& temperature_bin_exact_failures_used,
                           const std::vector<std::vector<size_t>>& temperature_bin_shadow_audits_used,
                           const std::vector<std::vector<size_t>>& temperature_bin_shadow_selected_decision_mismatches_used,
                           const std::vector<std::vector<size_t>>& temperature_bin_shadow_exact_failures_used,
                           const std::vector<std::vector<size_t>>& temperature_bin_shadow_selected_source_audits_used,
                           const std::vector<std::vector<size_t>>& temperature_bin_shadow_selected_source_exact_failures_used)
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
    require_size(primary_local_calibration_lookups_used.size(), "primary_local_calibration_lookups_used");
    require_size(primary_side_calibration_lookups_used.size(), "primary_side_calibration_lookups_used");
    require_size(primary_global_calibration_lookups_used.size(), "primary_global_calibration_lookups_used");
    require_size(primary_insufficient_calibration_lookups_used.size(), "primary_insufficient_calibration_lookups_used");
    require_size(aux_local_calibration_lookups_used.size(), "aux_local_calibration_lookups_used");
    require_size(aux_side_calibration_lookups_used.size(), "aux_side_calibration_lookups_used");
    require_size(aux_global_calibration_lookups_used.size(), "aux_global_calibration_lookups_used");
    require_size(aux_insufficient_calibration_lookups_used.size(), "aux_insufficient_calibration_lookups_used");
    require_size(mean_primary_ci_width_used.size(), "mean_primary_ci_width_used");
    require_size(mean_aux_ci_width_used.size(), "mean_aux_ci_width_used");
    require_size(mean_selected_ci_width_used.size(), "mean_selected_ci_width_used");
    require_size(mean_scanned_aux_levels_used.size(), "mean_scanned_aux_levels_used");
    require_size(shadow_audits_used.size(), "shadow_audits_used");
    require_size(shadow_exact_failures_used.size(), "shadow_exact_failures_used");
    require_size(shadow_primary_decision_mismatches_used.size(), "shadow_primary_decision_mismatches_used");
    require_size(shadow_selected_decision_mismatches_used.size(), "shadow_selected_decision_mismatches_used");
    require_size(shadow_selected_interval_hits_used.size(), "shadow_selected_interval_hits_used");
    require_size(shadow_full_scan_better_swaps_used.size(), "shadow_full_scan_better_swaps_used");
    require_size(mean_shadow_selected_regret_used.size(), "mean_shadow_selected_regret_used");
    require_size(mean_shadow_primary_abs_delta_error_used.size(), "mean_shadow_primary_abs_delta_error_used");
    require_size(mean_shadow_primary_delta_bias_used.size(), "mean_shadow_primary_delta_bias_used");
    require_size(mean_shadow_primary_p_bias_used.size(), "mean_shadow_primary_p_bias_used");
    require_size(mean_shadow_selected_abs_delta_error_used.size(), "mean_shadow_selected_abs_delta_error_used");
    require_size(mean_shadow_selected_abs_p_error_used.size(), "mean_shadow_selected_abs_p_error_used");
    require_size(mean_shadow_selected_delta_bias_used.size(), "mean_shadow_selected_delta_bias_used");
    require_size(mean_shadow_selected_p_bias_used.size(), "mean_shadow_selected_p_bias_used");
    require_size(shadow_selected_source_audits_used.size(), "shadow_selected_source_audits_used");
    require_size(shadow_selected_source_exact_failures_used.size(), "shadow_selected_source_exact_failures_used");
    require_size(shadow_selected_source_decision_mismatches_used.size(), "shadow_selected_source_decision_mismatches_used");
    require_size(shadow_selected_source_interval_hits_used.size(), "shadow_selected_source_interval_hits_used");
    require_size(shadow_selected_source_full_scan_better_swaps_used.size(), "shadow_selected_source_full_scan_better_swaps_used");
    require_size(mean_shadow_selected_source_regret_used.size(), "mean_shadow_selected_source_regret_used");
    require_size(mean_shadow_selected_source_abs_delta_error_used.size(), "mean_shadow_selected_source_abs_delta_error_used");
    require_size(mean_shadow_selected_source_abs_p_error_used.size(), "mean_shadow_selected_source_abs_p_error_used");
    require_size(mean_shadow_selected_source_delta_bias_used.size(), "mean_shadow_selected_source_delta_bias_used");
    require_size(mean_shadow_selected_source_p_bias_used.size(), "mean_shadow_selected_source_p_bias_used");
    require_size(temperature_bin_total_swaps_used.size(), "temperature_bin_total_swaps_used");
    require_size(temperature_bin_accepted_swaps_used.size(), "temperature_bin_accepted_swaps_used");
    require_size(temperature_bin_primary_resolved_swaps_used.size(), "temperature_bin_primary_resolved_swaps_used");
    require_size(temperature_bin_aux_resolved_swaps_used.size(), "temperature_bin_aux_resolved_swaps_used");
    require_size(temperature_bin_exact_evals_used.size(), "temperature_bin_exact_evals_used");
    require_size(temperature_bin_exact_unavailable_swaps_used.size(), "temperature_bin_exact_unavailable_swaps_used");
    require_size(temperature_bin_exact_failures_used.size(), "temperature_bin_exact_failures_used");
    require_size(temperature_bin_shadow_audits_used.size(), "temperature_bin_shadow_audits_used");
    require_size(temperature_bin_shadow_selected_decision_mismatches_used.size(), "temperature_bin_shadow_selected_decision_mismatches_used");
    require_size(temperature_bin_shadow_exact_failures_used.size(), "temperature_bin_shadow_exact_failures_used");
    require_size(temperature_bin_shadow_selected_source_audits_used.size(), "temperature_bin_shadow_selected_source_audits_used");
    require_size(temperature_bin_shadow_selected_source_exact_failures_used.size(), "temperature_bin_shadow_selected_source_exact_failures_used");

    for (size_t i = 0; i < n; ++i) {
        for (const auto* bins : {&temperature_bin_total_swaps_used, &temperature_bin_accepted_swaps_used,
                                 &temperature_bin_primary_resolved_swaps_used, &temperature_bin_aux_resolved_swaps_used,
                                 &temperature_bin_exact_evals_used, &temperature_bin_exact_unavailable_swaps_used,
                                 &temperature_bin_exact_failures_used, &temperature_bin_shadow_audits_used,
                                 &temperature_bin_shadow_selected_decision_mismatches_used,
                                         &temperature_bin_shadow_exact_failures_used}) {
            if (bins->at(i).size() != TEMPERATURE_BIN_COUNT) {
                throw std::logic_error("temperature-bin diagnostics have inconsistent width");
            }
        }

        for (const auto* source_temp_counts : {&temperature_bin_shadow_selected_source_audits_used,
                                              &temperature_bin_shadow_selected_source_exact_failures_used}) {
            if (source_temp_counts->at(i).size() != TEMPERATURE_BIN_COUNT * CALIBRATION_LOOKUP_SOURCE_COUNT) {
                throw std::logic_error("temperature-bin shadow source diagnostics have inconsistent width");
            }
        }

        for (const auto* source_counts : {&shadow_selected_source_audits_used,
                                          &shadow_selected_source_exact_failures_used,
                                          &shadow_selected_source_decision_mismatches_used,
                                          &shadow_selected_source_interval_hits_used,
                                          &shadow_selected_source_full_scan_better_swaps_used}) {
            if (source_counts->at(i).size() != CALIBRATION_LOOKUP_SOURCE_COUNT) {
                throw std::logic_error("shadow selected source diagnostics have inconsistent width");
            }
        }

        for (const auto* source_means : {&mean_shadow_selected_source_regret_used,
                                         &mean_shadow_selected_source_abs_delta_error_used,
                                         &mean_shadow_selected_source_abs_p_error_used,
                                         &mean_shadow_selected_source_delta_bias_used,
                                         &mean_shadow_selected_source_p_bias_used}) {
            if (source_means->at(i).size() != CALIBRATION_LOOKUP_SOURCE_COUNT) {
                throw std::logic_error("shadow selected source mean diagnostics have inconsistent width");
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
    size_t primary_local_calibration_lookups = 0;
    size_t primary_side_calibration_lookups = 0;
    size_t primary_global_calibration_lookups = 0;
    size_t primary_insufficient_calibration_lookups = 0;
    size_t aux_local_calibration_lookups = 0;
    size_t aux_side_calibration_lookups = 0;
    size_t aux_global_calibration_lookups = 0;
    size_t aux_insufficient_calibration_lookups = 0;
    double primary_ci_width_sum = 0.0;
    double aux_ci_width_sum = 0.0;
    double selected_ci_width_sum = 0.0;
    size_t primary_ci_width_observations = 0;
    size_t aux_ci_width_observations = 0;
    size_t selected_ci_width_observations = 0;
    double scanned_aux_levels_sum = 0.0;
    size_t scanned_aux_levels_observations = 0;
    size_t shadow_audits = 0;
    size_t shadow_exact_failures = 0;
    size_t shadow_primary_decision_mismatches = 0;
    size_t shadow_selected_decision_mismatches = 0;
    size_t shadow_selected_interval_hits = 0;
    size_t shadow_full_scan_better_swaps = 0;
    double shadow_selected_regret_sum = 0.0;
    double shadow_primary_abs_delta_error_sum = 0.0;
    double shadow_primary_delta_bias_sum = 0.0;
    double shadow_primary_p_bias_sum = 0.0;
    double shadow_selected_abs_delta_error_sum = 0.0;
    double shadow_selected_abs_p_error_sum = 0.0;
    double shadow_selected_delta_bias_sum = 0.0;
    double shadow_selected_p_bias_sum = 0.0;
    size_t shadow_metric_observations = 0;
    std::vector<size_t> shadow_selected_source_audits;
    std::vector<size_t> shadow_selected_source_exact_failures;
    std::vector<size_t> shadow_selected_source_decision_mismatches;
    std::vector<size_t> shadow_selected_source_interval_hits;
    std::vector<size_t> shadow_selected_source_full_scan_better_swaps;
    std::vector<double> shadow_selected_source_regret_sum;
    std::vector<double> shadow_selected_source_abs_delta_error_sum;
    std::vector<double> shadow_selected_source_abs_p_error_sum;
    std::vector<double> shadow_selected_source_delta_bias_sum;
    std::vector<double> shadow_selected_source_p_bias_sum;
    std::vector<size_t> shadow_selected_source_metric_observations;
    std::vector<size_t> temperature_bin_total_swaps;
    std::vector<size_t> temperature_bin_accepted_swaps;
    std::vector<size_t> temperature_bin_primary_resolved_swaps;
    std::vector<size_t> temperature_bin_aux_resolved_swaps;
    std::vector<size_t> temperature_bin_exact_evals;
    std::vector<size_t> temperature_bin_exact_unavailable_swaps;
    std::vector<size_t> temperature_bin_exact_failures;
    std::vector<size_t> temperature_bin_shadow_audits;
    std::vector<size_t> temperature_bin_shadow_selected_decision_mismatches;
    std::vector<size_t> temperature_bin_shadow_exact_failures;
    std::vector<size_t> temperature_bin_shadow_selected_source_audits;
    std::vector<size_t> temperature_bin_shadow_selected_source_exact_failures;
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

std::vector<double> filter_abs_errors(const std::vector<calibration_record>& records,
                                      const std::function<bool(const calibration_record&)>& predicate)
{
    std::vector<double> values;
    values.reserve(records.size());
    for (const auto& record : records) {
        if (predicate(record)) {
            values.push_back(record.abs_error);
        }
    }
    return values;
}

std::vector<double> filter_abs_deltas(const std::vector<calibration_record>& records,
                                      const std::function<bool(const calibration_record&)>& predicate)
{
    std::vector<double> values;
    values.reserve(records.size());
    for (const auto& record : records) {
        if (predicate(record)) {
            values.push_back(std::abs(record.cheap_delta));
        }
    }
    return values;
}

calibration_lookup_result adaptive_residual_bins::half_width(double cheap_delta, double calibration_quantile,
                                                             size_t min_calibration_samples) const
{
    calibration_lookup_result result;
    if (records.size() < min_calibration_samples) {
        return result;
    }

    const std::vector<double> abs_deltas = filter_abs_deltas(records, [](const calibration_record&) {
        return true;
    });
    const double center_abs_threshold = empirical_quantile(abs_deltas, ADAPTIVE_BIN_CENTER_ABS_QUANTILE);
    const bool in_center = std::abs(cheap_delta) <= center_abs_threshold;

    std::vector<double> local_errors;
    if (in_center) {
        local_errors = filter_abs_errors(records, [center_abs_threshold](const calibration_record& record) {
            return std::abs(record.cheap_delta) <= center_abs_threshold;
        });
    } else if (cheap_delta < 0.0) {
        const std::vector<double> negative_abs = filter_abs_deltas(records, [center_abs_threshold](const calibration_record& record) {
            return record.cheap_delta < -center_abs_threshold;
        });
        if (!negative_abs.empty()) {
            const double side_split = empirical_quantile(negative_abs, ADAPTIVE_BIN_SIDE_SPLIT_QUANTILE);
            const bool use_near = std::abs(cheap_delta) <= side_split;
            local_errors = filter_abs_errors(records, [center_abs_threshold, side_split, use_near](const calibration_record& record) {
                if (record.cheap_delta >= -center_abs_threshold) {
                    return false;
                }
                const double abs_delta = std::abs(record.cheap_delta);
                return use_near ? (abs_delta <= side_split) : (abs_delta > side_split);
            });
        }
    } else {
        const std::vector<double> positive_abs = filter_abs_deltas(records, [center_abs_threshold](const calibration_record& record) {
            return record.cheap_delta > center_abs_threshold;
        });
        if (!positive_abs.empty()) {
            const double side_split = empirical_quantile(positive_abs, ADAPTIVE_BIN_SIDE_SPLIT_QUANTILE);
            const bool use_near = std::abs(cheap_delta) <= side_split;
            local_errors = filter_abs_errors(records, [center_abs_threshold, side_split, use_near](const calibration_record& record) {
                if (record.cheap_delta <= center_abs_threshold) {
                    return false;
                }
                const double abs_delta = std::abs(record.cheap_delta);
                return use_near ? (abs_delta <= side_split) : (abs_delta > side_split);
            });
        }
    }
    if (local_errors.size() >= min_calibration_samples) {
        result.half_width = empirical_quantile(local_errors, calibration_quantile);
        result.sample_count = local_errors.size();
        result.source = calibration_lookup_source::local;
        return result;
    }

    std::vector<double> side_errors;
    if (cheap_delta < 0.0) {
        side_errors = filter_abs_errors(records, [](const calibration_record& record) {
            return record.cheap_delta < 0.0;
        });
    } else if (cheap_delta > 0.0) {
        side_errors = filter_abs_errors(records, [](const calibration_record& record) {
            return record.cheap_delta > 0.0;
        });
    }
    if (side_errors.size() >= min_calibration_samples) {
        result.half_width = empirical_quantile(side_errors, calibration_quantile);
        result.sample_count = side_errors.size();
        result.source = calibration_lookup_source::side;
        return result;
    }

    const std::vector<double> global_errors = filter_abs_errors(records, [](const calibration_record&) {
        return true;
    });
    if (global_errors.size() >= min_calibration_samples) {
        result.half_width = empirical_quantile(global_errors, calibration_quantile);
        result.sample_count = global_errors.size();
        result.source = calibration_lookup_source::global;
        return result;
    }

    return result;
}

double calibration_half_width(const std::vector<double>& residuals, double calibration_quantile,
                              size_t min_calibration_samples)
{
    if (residuals.size() < min_calibration_samples) {
        return std::numeric_limits<double>::infinity();
    }
    const double clamped_quantile = std::min(1.0, std::max(0.0, calibration_quantile));
    return empirical_quantile(residuals, clamped_quantile);
}

struct decision_interval {
    double lower_probability = 0.0;
    double upper_probability = 1.0;
    double width = 1.0;
    double point_probability = 0.0;
    bool resolved = false;
    bool accept = false;
};

decision_interval probability_interval_from_delta(double delta, double delta_half_width,
                                                  double temperature, double draw)
{
    decision_interval interval;
    interval.point_probability = acceptance_probability(delta, temperature);
    interval.accept = (draw < interval.point_probability);

    if (!std::isfinite(delta_half_width)) {
        return interval;
    }

    const double delta_lower = delta - delta_half_width;
    const double delta_upper = delta + delta_half_width;

    if (delta_upper <= 0.0) {
        interval.lower_probability = 1.0;
        interval.upper_probability = 1.0;
        interval.width = 0.0;
        interval.resolved = true;
        interval.accept = true;
        return interval;
    }

    interval.lower_probability = acceptance_probability(std::max(0.0, delta_upper), temperature);
    interval.upper_probability = (delta_lower <= 0.0)
        ? 1.0
        : acceptance_probability(delta_lower, temperature);
    interval.width = interval.upper_probability - interval.lower_probability;

    if (draw < interval.lower_probability) {
        interval.resolved = true;
        interval.accept = true;
    } else if (draw > interval.upper_probability) {
        interval.resolved = true;
        interval.accept = false;
    }

    return interval;
}

bool interval_contains_probability(const decision_interval& interval, double probability)
{
    constexpr double eps = 1e-12;
    return probability + eps >= interval.lower_probability &&
           probability <= interval.upper_probability + eps;
}

void push_bounded(std::vector<double>& values, double value, size_t max_calibration_history)
{
    if (values.size() >= max_calibration_history) {
        values.erase(values.begin());
    }
    values.push_back(value);
}

size_t calibration_lookup_source_index(calibration_lookup_source source)
{
    switch (source) {
    case calibration_lookup_source::local:
        return 0;
    case calibration_lookup_source::side:
        return 1;
    case calibration_lookup_source::global:
        return 2;
    case calibration_lookup_source::insufficient:
        return 3;
    }
    throw std::logic_error("unknown calibration lookup source");
}

void record_lookup_source(calibration_lookup_source source,
                          size_t& local_count,
                          size_t& side_count,
                          size_t& global_count,
                          size_t& insufficient_count)
{
    switch (source) {
    case calibration_lookup_source::local:
        ++local_count;
        break;
    case calibration_lookup_source::side:
        ++side_count;
        break;
    case calibration_lookup_source::global:
        ++global_count;
        break;
    case calibration_lookup_source::insufficient:
        ++insufficient_count;
        break;
    }
}

void record_lookup_source(calibration_lookup_source source,
                          std::vector<size_t>& counts)
{
    counts.at(calibration_lookup_source_index(source))++;
}

void record_temperature_source_lookup(calibration_lookup_source source,
                                      size_t temperature_bin,
                                      std::vector<size_t>& counts)
{
    const size_t index = temperature_bin * CALIBRATION_LOOKUP_SOURCE_COUNT +
                         calibration_lookup_source_index(source);
    counts.at(index)++;
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
    primary_local_calibration_lookups_used.clear();
    primary_side_calibration_lookups_used.clear();
    primary_global_calibration_lookups_used.clear();
    primary_insufficient_calibration_lookups_used.clear();
    aux_local_calibration_lookups_used.clear();
    aux_side_calibration_lookups_used.clear();
    aux_global_calibration_lookups_used.clear();
    aux_insufficient_calibration_lookups_used.clear();
    mean_primary_ci_width_used.clear();
    mean_aux_ci_width_used.clear();
    mean_selected_ci_width_used.clear();
    mean_scanned_aux_levels_used.clear();
    shadow_audits_used.clear();
    shadow_exact_failures_used.clear();
    shadow_primary_decision_mismatches_used.clear();
    shadow_selected_decision_mismatches_used.clear();
    shadow_selected_interval_hits_used.clear();
    shadow_full_scan_better_swaps_used.clear();
    mean_shadow_selected_regret_used.clear();
    mean_shadow_primary_abs_delta_error_used.clear();
    mean_shadow_primary_delta_bias_used.clear();
    mean_shadow_primary_p_bias_used.clear();
    mean_shadow_selected_abs_delta_error_used.clear();
    mean_shadow_selected_abs_p_error_used.clear();
    mean_shadow_selected_delta_bias_used.clear();
    mean_shadow_selected_p_bias_used.clear();
    shadow_selected_source_audits_used.clear();
    shadow_selected_source_exact_failures_used.clear();
    shadow_selected_source_decision_mismatches_used.clear();
    shadow_selected_source_interval_hits_used.clear();
    shadow_selected_source_full_scan_better_swaps_used.clear();
    mean_shadow_selected_source_regret_used.clear();
    mean_shadow_selected_source_abs_delta_error_used.clear();
    mean_shadow_selected_source_abs_p_error_used.clear();
    mean_shadow_selected_source_delta_bias_used.clear();
    mean_shadow_selected_source_p_bias_used.clear();
    temperature_bin_total_swaps_used.clear();
    temperature_bin_accepted_swaps_used.clear();
    temperature_bin_primary_resolved_swaps_used.clear();
    temperature_bin_aux_resolved_swaps_used.clear();
    temperature_bin_exact_evals_used.clear();
    temperature_bin_exact_unavailable_swaps_used.clear();
    temperature_bin_exact_failures_used.clear();
    temperature_bin_shadow_audits_used.clear();
    temperature_bin_shadow_selected_decision_mismatches_used.clear();
    temperature_bin_shadow_exact_failures_used.clear();
    temperature_bin_shadow_selected_source_audits_used.clear();
    temperature_bin_shadow_selected_source_exact_failures_used.clear();

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
                adaptive_residual_bins primary_delta_residuals(max_calibration_history);
                std::vector<adaptive_residual_bins> aux_delta_residuals;
                aux_delta_residuals.reserve(local_test->aux_statistic_levels());
                for (size_t level = 0; level < local_test->aux_statistic_levels(); ++level) {
                    aux_delta_residuals.emplace_back(max_calibration_history);
                }
                result.aux_level_resolved_swaps.assign(local_test->aux_statistic_levels(), 0);
                result.temperature_bin_total_swaps.assign(TEMPERATURE_BIN_COUNT, 0);
                result.temperature_bin_accepted_swaps.assign(TEMPERATURE_BIN_COUNT, 0);
                result.temperature_bin_primary_resolved_swaps.assign(TEMPERATURE_BIN_COUNT, 0);
                result.temperature_bin_aux_resolved_swaps.assign(TEMPERATURE_BIN_COUNT, 0);
                result.temperature_bin_exact_evals.assign(TEMPERATURE_BIN_COUNT, 0);
                result.temperature_bin_exact_unavailable_swaps.assign(TEMPERATURE_BIN_COUNT, 0);
                result.temperature_bin_exact_failures.assign(TEMPERATURE_BIN_COUNT, 0);
                result.temperature_bin_shadow_audits.assign(TEMPERATURE_BIN_COUNT, 0);
                result.temperature_bin_shadow_selected_decision_mismatches.assign(TEMPERATURE_BIN_COUNT, 0);
                result.temperature_bin_shadow_exact_failures.assign(TEMPERATURE_BIN_COUNT, 0);
                result.temperature_bin_shadow_selected_source_audits.assign(
                    TEMPERATURE_BIN_COUNT * CALIBRATION_LOOKUP_SOURCE_COUNT,
                    0
                );
                result.temperature_bin_shadow_selected_source_exact_failures.assign(
                    TEMPERATURE_BIN_COUNT * CALIBRATION_LOOKUP_SOURCE_COUNT,
                    0
                );
                result.shadow_selected_source_audits.assign(CALIBRATION_LOOKUP_SOURCE_COUNT, 0);
                result.shadow_selected_source_exact_failures.assign(CALIBRATION_LOOKUP_SOURCE_COUNT, 0);
                result.shadow_selected_source_decision_mismatches.assign(CALIBRATION_LOOKUP_SOURCE_COUNT, 0);
                result.shadow_selected_source_interval_hits.assign(CALIBRATION_LOOKUP_SOURCE_COUNT, 0);
                result.shadow_selected_source_full_scan_better_swaps.assign(CALIBRATION_LOOKUP_SOURCE_COUNT, 0);
                result.shadow_selected_source_regret_sum.assign(CALIBRATION_LOOKUP_SOURCE_COUNT, 0.0);
                result.shadow_selected_source_abs_delta_error_sum.assign(CALIBRATION_LOOKUP_SOURCE_COUNT, 0.0);
                result.shadow_selected_source_abs_p_error_sum.assign(CALIBRATION_LOOKUP_SOURCE_COUNT, 0.0);
                result.shadow_selected_source_delta_bias_sum.assign(CALIBRATION_LOOKUP_SOURCE_COUNT, 0.0);
                result.shadow_selected_source_p_bias_sum.assign(CALIBRATION_LOOKUP_SOURCE_COUNT, 0.0);
                result.shadow_selected_source_metric_observations.assign(CALIBRATION_LOOKUP_SOURCE_COUNT, 0);
                result.total_swaps = iterations;
                while (local_test->subsample_size() < curr_size) {
                    local_test->add_one();
                }
                std::uniform_real_distribution<double> random_unif(0.0, 1.0);
                double score = local_test->get_normality_statistic();
                for (size_t k = 0; k < iterations; k++) {
                    t = c * t;
                    const size_t temp_bin = temperature_bin_index(k, iterations);
                    result.temperature_bin_total_swaps[temp_bin]++;
                    local_test->swap_once();
                    double new_score = local_test->get_normality_statistic();

                    const double delta_primary = new_score - score;
                    const bool primary_improving = delta_primary <= 0.0;
                    const double p_primary = acceptance_probability(delta_primary, t);
                    const double acceptance_draw = random_unif(mersenne_wheel);

                    bool accept = true;

                    if (local_test->has_aux_statistic()) {
                        const bool calibration_ready = primary_delta_residuals.size() >= min_calibration_samples;
                        if (!calibration_ready) {
                            result.primary_ci_width_sum += 1.0;
                            result.primary_ci_width_observations++;
                            result.aux_ci_width_sum += 1.0;
                            result.aux_ci_width_observations++;
                            result.selected_ci_width_sum += 1.0;
                            result.selected_ci_width_observations++;
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

                                    const std::vector<double> aux_deltas = local_test->aux_deltas_last_swap();
                                    const double delta_exact = local_test->exact_delta_last_swap();
                                    const double p_exact = acceptance_probability(delta_exact, t);
                                    primary_delta_residuals.observe(delta_primary, delta_exact);
                                    for (size_t level = 0; level < aux_deltas.size(); ++level) {
                                        aux_delta_residuals[level].observe(aux_deltas[level], delta_exact);
                                    }
                                    accept = (acceptance_draw < p_exact);
                                } catch (...) {
                                    result.exact_failures++;
                                    result.temperature_bin_exact_failures[temp_bin]++;
                                    accept = (acceptance_draw < p_primary);
                                }
                            } else {
                                result.exact_unavailable_swaps++;
                                result.temperature_bin_exact_unavailable_swaps[temp_bin]++;
                                accept = (acceptance_draw < p_primary);
                            }
                        } else {
                            const calibration_lookup_result primary_lookup = primary_delta_residuals.half_width(
                                delta_primary, calibration_quantile, min_calibration_samples);
                            record_lookup_source(primary_lookup.source,
                                                 result.primary_local_calibration_lookups,
                                                 result.primary_side_calibration_lookups,
                                                 result.primary_global_calibration_lookups,
                                                 result.primary_insufficient_calibration_lookups);
                            const decision_interval primary_interval = probability_interval_from_delta(
                                delta_primary,
                                primary_lookup.half_width,
                                t,
                                acceptance_draw);

                            double representative_aux_width = 1.0;
                            double selected_width = primary_interval.width;
                            bool cheap_accept = primary_interval.accept;
                            double fallback_width = primary_interval.width;
                            bool fallback_accept = primary_interval.accept;
                            bool resolved_by_ci = false;
                            bool resolved_by_primary = false;
                            bool resolved_by_aux = false;
                            bool observed_aux_width = false;
                            size_t scanned_aux_levels = 0;
                            size_t resolved_aux_level = local_test->aux_statistic_levels();
                            double selected_delta = delta_primary;
                            decision_interval selected_interval = primary_interval;
                            calibration_lookup_source selected_lookup_source = primary_lookup.source;

                            if (primary_interval.resolved && primary_interval.width <= ci_width_threshold) {
                                resolved_by_ci = true;
                                resolved_by_primary = true;
                                cheap_accept = primary_interval.accept;
                                selected_width = primary_interval.width;
                                selected_lookup_source = primary_lookup.source;
                            } else {
                                local_test->scan_aux_deltas_last_swap([&](size_t level, double delta_aux) {
                                    scanned_aux_levels++;
                                    const calibration_lookup_result aux_lookup = aux_delta_residuals[level].half_width(
                                        delta_aux, calibration_quantile, min_calibration_samples);
                                    record_lookup_source(aux_lookup.source,
                                                         result.aux_local_calibration_lookups,
                                                         result.aux_side_calibration_lookups,
                                                         result.aux_global_calibration_lookups,
                                                         result.aux_insufficient_calibration_lookups);
                                    const decision_interval aux_interval = probability_interval_from_delta(
                                        delta_aux,
                                        aux_lookup.half_width,
                                        t,
                                        acceptance_draw);
                                    observed_aux_width = true;
                                    representative_aux_width = std::min(representative_aux_width, aux_interval.width);
                                    if (aux_interval.width < fallback_width) {
                                        fallback_width = aux_interval.width;
                                        fallback_accept = aux_interval.accept;
                                    }
                                    if (aux_interval.resolved && aux_interval.width <= ci_width_threshold) {
                                        resolved_by_ci = true;
                                        resolved_by_aux = true;
                                        resolved_aux_level = level;
                                        cheap_accept = aux_interval.accept;
                                        selected_width = aux_interval.width;
                                        selected_delta = delta_aux;
                                        selected_interval = aux_interval;
                                        selected_lookup_source = aux_lookup.source;
                                        return false;
                                    }
                                    return true;
                                });
                            }

                            result.primary_ci_width_sum += primary_interval.width;
                            result.primary_ci_width_observations++;
                            result.scanned_aux_levels_sum += static_cast<double>(scanned_aux_levels);
                            result.scanned_aux_levels_observations++;
                            if (observed_aux_width) {
                                result.aux_ci_width_sum += representative_aux_width;
                                result.aux_ci_width_observations++;
                            }

                            if (!resolved_by_ci) {
                                selected_width = fallback_width;
                                cheap_accept = fallback_accept;
                            }
                            result.selected_ci_width_sum += selected_width;
                            result.selected_ci_width_observations++;

                            const bool should_audit = resolved_by_ci && local_test->last_swap_has_equal_effect_size() &&
                                (primary_delta_residuals.size() < min_calibration_samples ||
                                 random_unif(mersenne_wheel) < HYBRID_RESOLVED_AUDIT_RATE);
                            const bool need_exact = !resolved_by_ci || should_audit;

                            if (need_exact) {
                                if (!resolved_by_ci) {
                                    result.uncertain_swaps++;
                                }
                                if (local_test->last_swap_has_equal_effect_size()) {
                                    try {
                                        result.exact_evals++;
                                        result.temperature_bin_exact_evals[temp_bin]++;
                                        if (primary_improving) {
                                            result.exact_evals_on_improving++;
                                        } else {
                                            result.exact_evals_on_worsening++;
                                        }

                                        const std::vector<double> aux_deltas = local_test->aux_deltas_last_swap();
                                        const double delta_exact = local_test->exact_delta_last_swap();
                                        const double p_exact = acceptance_probability(delta_exact, t);
                                        primary_delta_residuals.observe(delta_primary, delta_exact);
                                        for (size_t level = 0; level < aux_deltas.size(); ++level) {
                                            aux_delta_residuals[level].observe(aux_deltas[level], delta_exact);
                                        }
                                        accept = (acceptance_draw < p_exact);
                                    } catch (...) {
                                        result.exact_failures++;
                                        result.temperature_bin_exact_failures[temp_bin]++;

                                         if (resolved_by_ci) {
                                            result.ci_resolved_swaps++;
                                            if (resolved_by_primary) {
                                                result.primary_resolved_swaps++;
                                                result.temperature_bin_primary_resolved_swaps[temp_bin]++;
                                            }
                                            if (resolved_by_aux && resolved_aux_level < result.aux_level_resolved_swaps.size()) {
                                                result.aux_level_resolved_swaps[resolved_aux_level]++;
                                                result.temperature_bin_aux_resolved_swaps[temp_bin]++;
                                            }
                                        }
                                        accept = cheap_accept;
                                    }
                                } else {
                                    result.exact_unavailable_swaps++;
                                    result.temperature_bin_exact_unavailable_swaps[temp_bin]++;
                                    accept = cheap_accept;
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
                                accept = cheap_accept;

                                const bool should_shadow_audit = local_test->last_swap_has_equal_effect_size() &&
                                    random_unif(mersenne_wheel) < HYBRID_SHADOW_AUDIT_RATE;
                                if (should_shadow_audit) {
                                    result.shadow_audits++;
                                    result.temperature_bin_shadow_audits[temp_bin]++;
                                    record_lookup_source(selected_lookup_source, result.shadow_selected_source_audits);
                                    record_temperature_source_lookup(
                                        selected_lookup_source,
                                        temp_bin,
                                        result.temperature_bin_shadow_selected_source_audits
                                    );
                                    const size_t selected_source_index = calibration_lookup_source_index(selected_lookup_source);
                                    try {
                                        const std::vector<double> aux_deltas = local_test->aux_deltas_last_swap();
                                        const double delta_exact = local_test->exact_delta_last_swap();
                                        const double p_exact = acceptance_probability(delta_exact, t);
                                        const bool exact_accept = (acceptance_draw < p_exact);

                                        result.shadow_metric_observations++;
                                        result.shadow_primary_abs_delta_error_sum += std::abs(delta_exact - delta_primary);
                                        result.shadow_primary_delta_bias_sum += (delta_primary - delta_exact);
                                        result.shadow_primary_p_bias_sum += (acceptance_probability(delta_primary, t) - p_exact);
                                        result.shadow_selected_abs_delta_error_sum += std::abs(delta_exact - selected_delta);
                                        result.shadow_selected_abs_p_error_sum += std::abs(p_exact - acceptance_probability(selected_delta, t));
                                        result.shadow_selected_delta_bias_sum += (selected_delta - delta_exact);
                                        result.shadow_selected_p_bias_sum += (acceptance_probability(selected_delta, t) - p_exact);
                                        result.shadow_selected_source_metric_observations[selected_source_index]++;
                                        result.shadow_selected_source_abs_delta_error_sum[selected_source_index] +=
                                            std::abs(delta_exact - selected_delta);
                                        result.shadow_selected_source_abs_p_error_sum[selected_source_index] +=
                                            std::abs(p_exact - acceptance_probability(selected_delta, t));
                                        result.shadow_selected_source_delta_bias_sum[selected_source_index] +=
                                            (selected_delta - delta_exact);
                                        result.shadow_selected_source_p_bias_sum[selected_source_index] +=
                                            (acceptance_probability(selected_delta, t) - p_exact);

                                        if (primary_interval.accept != exact_accept) {
                                            result.shadow_primary_decision_mismatches++;
                                        }
                                        if (cheap_accept != exact_accept) {
                                            result.shadow_selected_decision_mismatches++;
                                            result.temperature_bin_shadow_selected_decision_mismatches[temp_bin]++;
                                            result.shadow_selected_source_decision_mismatches[selected_source_index]++;
                                        }
                                        if (interval_contains_probability(selected_interval, p_exact)) {
                                            result.shadow_selected_interval_hits++;
                                            result.shadow_selected_source_interval_hits[selected_source_index]++;
                                        }

                                        double full_best_width = std::numeric_limits<double>::infinity();
                                        if (primary_interval.resolved && primary_interval.width <= ci_width_threshold) {
                                            full_best_width = primary_interval.width;
                                        }
                                        for (size_t level = 0; level < aux_deltas.size(); ++level) {
                                            const calibration_lookup_result aux_lookup = aux_delta_residuals[level].half_width(
                                                aux_deltas[level], calibration_quantile, min_calibration_samples);
                                            const decision_interval aux_interval = probability_interval_from_delta(
                                                aux_deltas[level],
                                                aux_lookup.half_width,
                                                t,
                                                acceptance_draw);
                                            if (aux_interval.resolved && aux_interval.width <= ci_width_threshold) {
                                                full_best_width = std::min(full_best_width, aux_interval.width);
                                            }
                                        }

                                        const double regret = std::isfinite(full_best_width)
                                            ? std::max(0.0, selected_width - full_best_width)
                                            : 0.0;
                                        result.shadow_selected_regret_sum += regret;
                                        result.shadow_selected_source_regret_sum[selected_source_index] += regret;
                                        if (regret > 0.0) {
                                            result.shadow_full_scan_better_swaps++;
                                            result.shadow_selected_source_full_scan_better_swaps[selected_source_index]++;
                                        }
                                    } catch (...) {
                                        result.shadow_exact_failures++;
                                        result.temperature_bin_shadow_exact_failures[temp_bin]++;
                                        result.shadow_selected_source_exact_failures[selected_source_index]++;
                                        record_temperature_source_lookup(
                                            selected_lookup_source,
                                            temp_bin,
                                            result.temperature_bin_shadow_selected_source_exact_failures
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        accept = (acceptance_draw < p_primary);
                    }

                    if (!accept) {
                        local_test->swap_once(true);
                    } else {
                        result.temperature_bin_accepted_swaps[temp_bin]++;
                        score = new_score;
                    }
                }
                result.primary_calibration_points = primary_delta_residuals.size();
                result.aux_calibration_points = aux_delta_residuals.empty() ? 0 : aux_delta_residuals.back().size();
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
        size_t primary_local_calibration_lookups_curr = 0;
        size_t primary_side_calibration_lookups_curr = 0;
        size_t primary_global_calibration_lookups_curr = 0;
        size_t primary_insufficient_calibration_lookups_curr = 0;
        size_t aux_local_calibration_lookups_curr = 0;
        size_t aux_side_calibration_lookups_curr = 0;
        size_t aux_global_calibration_lookups_curr = 0;
        size_t aux_insufficient_calibration_lookups_curr = 0;
        double primary_ci_width_sum_curr = 0.0;
        double aux_ci_width_sum_curr = 0.0;
        double selected_ci_width_sum_curr = 0.0;
        size_t primary_ci_width_observations_curr = 0;
        size_t aux_ci_width_observations_curr = 0;
        size_t selected_ci_width_observations_curr = 0;
        double scanned_aux_levels_sum_curr = 0.0;
        size_t scanned_aux_levels_observations_curr = 0;
        size_t shadow_audits_curr = 0;
        size_t shadow_exact_failures_curr = 0;
        size_t shadow_primary_decision_mismatches_curr = 0;
        size_t shadow_selected_decision_mismatches_curr = 0;
        size_t shadow_selected_interval_hits_curr = 0;
        size_t shadow_full_scan_better_swaps_curr = 0;
        double shadow_selected_regret_sum_curr = 0.0;
        double shadow_primary_abs_delta_error_sum_curr = 0.0;
        double shadow_primary_delta_bias_sum_curr = 0.0;
        double shadow_primary_p_bias_sum_curr = 0.0;
        double shadow_selected_abs_delta_error_sum_curr = 0.0;
        double shadow_selected_abs_p_error_sum_curr = 0.0;
        double shadow_selected_delta_bias_sum_curr = 0.0;
        double shadow_selected_p_bias_sum_curr = 0.0;
        size_t shadow_metric_observations_curr = 0;
        std::vector<size_t> shadow_selected_source_audits_curr(CALIBRATION_LOOKUP_SOURCE_COUNT, 0);
        std::vector<size_t> shadow_selected_source_exact_failures_curr(CALIBRATION_LOOKUP_SOURCE_COUNT, 0);
        std::vector<size_t> shadow_selected_source_decision_mismatches_curr(CALIBRATION_LOOKUP_SOURCE_COUNT, 0);
        std::vector<size_t> shadow_selected_source_interval_hits_curr(CALIBRATION_LOOKUP_SOURCE_COUNT, 0);
        std::vector<size_t> shadow_selected_source_full_scan_better_swaps_curr(CALIBRATION_LOOKUP_SOURCE_COUNT, 0);
        std::vector<double> shadow_selected_source_regret_sum_curr(CALIBRATION_LOOKUP_SOURCE_COUNT, 0.0);
        std::vector<double> shadow_selected_source_abs_delta_error_sum_curr(CALIBRATION_LOOKUP_SOURCE_COUNT, 0.0);
        std::vector<double> shadow_selected_source_abs_p_error_sum_curr(CALIBRATION_LOOKUP_SOURCE_COUNT, 0.0);
        std::vector<double> shadow_selected_source_delta_bias_sum_curr(CALIBRATION_LOOKUP_SOURCE_COUNT, 0.0);
        std::vector<double> shadow_selected_source_p_bias_sum_curr(CALIBRATION_LOOKUP_SOURCE_COUNT, 0.0);
        std::vector<size_t> shadow_selected_source_metric_observations_curr(CALIBRATION_LOOKUP_SOURCE_COUNT, 0);
        std::vector<size_t> temperature_bin_total_swaps_curr(TEMPERATURE_BIN_COUNT, 0);
        std::vector<size_t> temperature_bin_accepted_swaps_curr(TEMPERATURE_BIN_COUNT, 0);
        std::vector<size_t> temperature_bin_primary_resolved_swaps_curr(TEMPERATURE_BIN_COUNT, 0);
        std::vector<size_t> temperature_bin_aux_resolved_swaps_curr(TEMPERATURE_BIN_COUNT, 0);
        std::vector<size_t> temperature_bin_exact_evals_curr(TEMPERATURE_BIN_COUNT, 0);
        std::vector<size_t> temperature_bin_exact_unavailable_swaps_curr(TEMPERATURE_BIN_COUNT, 0);
        std::vector<size_t> temperature_bin_exact_failures_curr(TEMPERATURE_BIN_COUNT, 0);
        std::vector<size_t> temperature_bin_shadow_audits_curr(TEMPERATURE_BIN_COUNT, 0);
        std::vector<size_t> temperature_bin_shadow_selected_decision_mismatches_curr(TEMPERATURE_BIN_COUNT, 0);
        std::vector<size_t> temperature_bin_shadow_exact_failures_curr(TEMPERATURE_BIN_COUNT, 0);
        std::vector<size_t> temperature_bin_shadow_selected_source_audits_curr(
            TEMPERATURE_BIN_COUNT * CALIBRATION_LOOKUP_SOURCE_COUNT,
            0
        );
        std::vector<size_t> temperature_bin_shadow_selected_source_exact_failures_curr(
            TEMPERATURE_BIN_COUNT * CALIBRATION_LOOKUP_SOURCE_COUNT,
            0
        );

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
            primary_local_calibration_lookups_curr += run.primary_local_calibration_lookups;
            primary_side_calibration_lookups_curr += run.primary_side_calibration_lookups;
            primary_global_calibration_lookups_curr += run.primary_global_calibration_lookups;
            primary_insufficient_calibration_lookups_curr += run.primary_insufficient_calibration_lookups;
            aux_local_calibration_lookups_curr += run.aux_local_calibration_lookups;
            aux_side_calibration_lookups_curr += run.aux_side_calibration_lookups;
            aux_global_calibration_lookups_curr += run.aux_global_calibration_lookups;
            aux_insufficient_calibration_lookups_curr += run.aux_insufficient_calibration_lookups;
            primary_ci_width_sum_curr += run.primary_ci_width_sum;
            aux_ci_width_sum_curr += run.aux_ci_width_sum;
            selected_ci_width_sum_curr += run.selected_ci_width_sum;
            primary_ci_width_observations_curr += run.primary_ci_width_observations;
            aux_ci_width_observations_curr += run.aux_ci_width_observations;
            selected_ci_width_observations_curr += run.selected_ci_width_observations;
            scanned_aux_levels_sum_curr += run.scanned_aux_levels_sum;
            scanned_aux_levels_observations_curr += run.scanned_aux_levels_observations;
            shadow_audits_curr += run.shadow_audits;
            shadow_exact_failures_curr += run.shadow_exact_failures;
            shadow_primary_decision_mismatches_curr += run.shadow_primary_decision_mismatches;
            shadow_selected_decision_mismatches_curr += run.shadow_selected_decision_mismatches;
            shadow_selected_interval_hits_curr += run.shadow_selected_interval_hits;
            shadow_full_scan_better_swaps_curr += run.shadow_full_scan_better_swaps;
            shadow_selected_regret_sum_curr += run.shadow_selected_regret_sum;
            shadow_primary_abs_delta_error_sum_curr += run.shadow_primary_abs_delta_error_sum;
            shadow_primary_delta_bias_sum_curr += run.shadow_primary_delta_bias_sum;
            shadow_primary_p_bias_sum_curr += run.shadow_primary_p_bias_sum;
            shadow_selected_abs_delta_error_sum_curr += run.shadow_selected_abs_delta_error_sum;
            shadow_selected_abs_p_error_sum_curr += run.shadow_selected_abs_p_error_sum;
            shadow_selected_delta_bias_sum_curr += run.shadow_selected_delta_bias_sum;
            shadow_selected_p_bias_sum_curr += run.shadow_selected_p_bias_sum;
            shadow_metric_observations_curr += run.shadow_metric_observations;
            for (size_t source = 0; source < CALIBRATION_LOOKUP_SOURCE_COUNT; ++source) {
                shadow_selected_source_audits_curr[source] += run.shadow_selected_source_audits[source];
                shadow_selected_source_exact_failures_curr[source] += run.shadow_selected_source_exact_failures[source];
                shadow_selected_source_decision_mismatches_curr[source] += run.shadow_selected_source_decision_mismatches[source];
                shadow_selected_source_interval_hits_curr[source] += run.shadow_selected_source_interval_hits[source];
                shadow_selected_source_full_scan_better_swaps_curr[source] += run.shadow_selected_source_full_scan_better_swaps[source];
                shadow_selected_source_regret_sum_curr[source] += run.shadow_selected_source_regret_sum[source];
                shadow_selected_source_abs_delta_error_sum_curr[source] += run.shadow_selected_source_abs_delta_error_sum[source];
                shadow_selected_source_abs_p_error_sum_curr[source] += run.shadow_selected_source_abs_p_error_sum[source];
                shadow_selected_source_delta_bias_sum_curr[source] += run.shadow_selected_source_delta_bias_sum[source];
                shadow_selected_source_p_bias_sum_curr[source] += run.shadow_selected_source_p_bias_sum[source];
                shadow_selected_source_metric_observations_curr[source] += run.shadow_selected_source_metric_observations[source];
            }
            for (size_t bin = 0; bin < TEMPERATURE_BIN_COUNT; ++bin) {
                temperature_bin_total_swaps_curr[bin] += run.temperature_bin_total_swaps[bin];
                temperature_bin_accepted_swaps_curr[bin] += run.temperature_bin_accepted_swaps[bin];
                temperature_bin_primary_resolved_swaps_curr[bin] += run.temperature_bin_primary_resolved_swaps[bin];
                temperature_bin_aux_resolved_swaps_curr[bin] += run.temperature_bin_aux_resolved_swaps[bin];
                temperature_bin_exact_evals_curr[bin] += run.temperature_bin_exact_evals[bin];
                temperature_bin_exact_unavailable_swaps_curr[bin] += run.temperature_bin_exact_unavailable_swaps[bin];
                temperature_bin_exact_failures_curr[bin] += run.temperature_bin_exact_failures[bin];
                temperature_bin_shadow_audits_curr[bin] += run.temperature_bin_shadow_audits[bin];
                temperature_bin_shadow_selected_decision_mismatches_curr[bin] += run.temperature_bin_shadow_selected_decision_mismatches[bin];
                temperature_bin_shadow_exact_failures_curr[bin] += run.temperature_bin_shadow_exact_failures[bin];
            }
            for (size_t index = 0; index < temperature_bin_shadow_selected_source_audits_curr.size(); ++index) {
                temperature_bin_shadow_selected_source_audits_curr[index] += run.temperature_bin_shadow_selected_source_audits[index];
                temperature_bin_shadow_selected_source_exact_failures_curr[index] += run.temperature_bin_shadow_selected_source_exact_failures[index];
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
        primary_local_calibration_lookups_used.push_back(primary_local_calibration_lookups_curr);
        primary_side_calibration_lookups_used.push_back(primary_side_calibration_lookups_curr);
        primary_global_calibration_lookups_used.push_back(primary_global_calibration_lookups_curr);
        primary_insufficient_calibration_lookups_used.push_back(primary_insufficient_calibration_lookups_curr);
        aux_local_calibration_lookups_used.push_back(aux_local_calibration_lookups_curr);
        aux_side_calibration_lookups_used.push_back(aux_side_calibration_lookups_curr);
        aux_global_calibration_lookups_used.push_back(aux_global_calibration_lookups_curr);
        aux_insufficient_calibration_lookups_used.push_back(aux_insufficient_calibration_lookups_curr);
        mean_primary_ci_width_used.push_back(
            (primary_ci_width_observations_curr == 0)
                ? NA_REAL
                : primary_ci_width_sum_curr / static_cast<double>(primary_ci_width_observations_curr));
        mean_aux_ci_width_used.push_back(
            (aux_ci_width_observations_curr == 0)
                ? NA_REAL
                : aux_ci_width_sum_curr / static_cast<double>(aux_ci_width_observations_curr));
        mean_selected_ci_width_used.push_back(
            (selected_ci_width_observations_curr == 0)
                ? NA_REAL
                : selected_ci_width_sum_curr / static_cast<double>(selected_ci_width_observations_curr));
        mean_scanned_aux_levels_used.push_back(
            (scanned_aux_levels_observations_curr == 0)
                ? NA_REAL
                : scanned_aux_levels_sum_curr / static_cast<double>(scanned_aux_levels_observations_curr));
        shadow_audits_used.push_back(shadow_audits_curr);
        shadow_exact_failures_used.push_back(shadow_exact_failures_curr);
        shadow_primary_decision_mismatches_used.push_back(shadow_primary_decision_mismatches_curr);
        shadow_selected_decision_mismatches_used.push_back(shadow_selected_decision_mismatches_curr);
        shadow_selected_interval_hits_used.push_back(shadow_selected_interval_hits_curr);
        shadow_full_scan_better_swaps_used.push_back(shadow_full_scan_better_swaps_curr);
        mean_shadow_selected_regret_used.push_back(
            (shadow_metric_observations_curr == 0)
                ? NA_REAL
                : shadow_selected_regret_sum_curr / static_cast<double>(shadow_metric_observations_curr));
        mean_shadow_primary_abs_delta_error_used.push_back(
            (shadow_metric_observations_curr == 0)
                ? NA_REAL
                : shadow_primary_abs_delta_error_sum_curr / static_cast<double>(shadow_metric_observations_curr));
        mean_shadow_primary_delta_bias_used.push_back(
            (shadow_metric_observations_curr == 0)
                ? NA_REAL
                : shadow_primary_delta_bias_sum_curr / static_cast<double>(shadow_metric_observations_curr));
        mean_shadow_primary_p_bias_used.push_back(
            (shadow_metric_observations_curr == 0)
                ? NA_REAL
                : shadow_primary_p_bias_sum_curr / static_cast<double>(shadow_metric_observations_curr));
        mean_shadow_selected_abs_delta_error_used.push_back(
            (shadow_metric_observations_curr == 0)
                ? NA_REAL
                : shadow_selected_abs_delta_error_sum_curr / static_cast<double>(shadow_metric_observations_curr));
        mean_shadow_selected_abs_p_error_used.push_back(
            (shadow_metric_observations_curr == 0)
                ? NA_REAL
                : shadow_selected_abs_p_error_sum_curr / static_cast<double>(shadow_metric_observations_curr));
        mean_shadow_selected_delta_bias_used.push_back(
            (shadow_metric_observations_curr == 0)
                ? NA_REAL
                : shadow_selected_delta_bias_sum_curr / static_cast<double>(shadow_metric_observations_curr));
        mean_shadow_selected_p_bias_used.push_back(
            (shadow_metric_observations_curr == 0)
                ? NA_REAL
                : shadow_selected_p_bias_sum_curr / static_cast<double>(shadow_metric_observations_curr));
        shadow_selected_source_audits_used.push_back(shadow_selected_source_audits_curr);
        shadow_selected_source_exact_failures_used.push_back(shadow_selected_source_exact_failures_curr);
        shadow_selected_source_decision_mismatches_used.push_back(shadow_selected_source_decision_mismatches_curr);
        shadow_selected_source_interval_hits_used.push_back(shadow_selected_source_interval_hits_curr);
        shadow_selected_source_full_scan_better_swaps_used.push_back(shadow_selected_source_full_scan_better_swaps_curr);
        std::vector<double> mean_shadow_selected_source_regret_curr(CALIBRATION_LOOKUP_SOURCE_COUNT, NA_REAL);
        std::vector<double> mean_shadow_selected_source_abs_delta_error_curr(CALIBRATION_LOOKUP_SOURCE_COUNT, NA_REAL);
        std::vector<double> mean_shadow_selected_source_abs_p_error_curr(CALIBRATION_LOOKUP_SOURCE_COUNT, NA_REAL);
        std::vector<double> mean_shadow_selected_source_delta_bias_curr(CALIBRATION_LOOKUP_SOURCE_COUNT, NA_REAL);
        std::vector<double> mean_shadow_selected_source_p_bias_curr(CALIBRATION_LOOKUP_SOURCE_COUNT, NA_REAL);
        for (size_t source = 0; source < CALIBRATION_LOOKUP_SOURCE_COUNT; ++source) {
            if (shadow_selected_source_metric_observations_curr[source] == 0) {
                continue;
            }
            const double denom = static_cast<double>(shadow_selected_source_metric_observations_curr[source]);
            mean_shadow_selected_source_regret_curr[source] =
                shadow_selected_source_regret_sum_curr[source] / denom;
            mean_shadow_selected_source_abs_delta_error_curr[source] =
                shadow_selected_source_abs_delta_error_sum_curr[source] / denom;
            mean_shadow_selected_source_abs_p_error_curr[source] =
                shadow_selected_source_abs_p_error_sum_curr[source] / denom;
            mean_shadow_selected_source_delta_bias_curr[source] =
                shadow_selected_source_delta_bias_sum_curr[source] / denom;
            mean_shadow_selected_source_p_bias_curr[source] =
                shadow_selected_source_p_bias_sum_curr[source] / denom;
        }
        mean_shadow_selected_source_regret_used.push_back(mean_shadow_selected_source_regret_curr);
        mean_shadow_selected_source_abs_delta_error_used.push_back(mean_shadow_selected_source_abs_delta_error_curr);
        mean_shadow_selected_source_abs_p_error_used.push_back(mean_shadow_selected_source_abs_p_error_curr);
        mean_shadow_selected_source_delta_bias_used.push_back(mean_shadow_selected_source_delta_bias_curr);
        mean_shadow_selected_source_p_bias_used.push_back(mean_shadow_selected_source_p_bias_curr);
        temperature_bin_total_swaps_used.push_back(temperature_bin_total_swaps_curr);
        temperature_bin_accepted_swaps_used.push_back(temperature_bin_accepted_swaps_curr);
        temperature_bin_primary_resolved_swaps_used.push_back(temperature_bin_primary_resolved_swaps_curr);
        temperature_bin_aux_resolved_swaps_used.push_back(temperature_bin_aux_resolved_swaps_curr);
        temperature_bin_exact_evals_used.push_back(temperature_bin_exact_evals_curr);
        temperature_bin_exact_unavailable_swaps_used.push_back(temperature_bin_exact_unavailable_swaps_curr);
        temperature_bin_exact_failures_used.push_back(temperature_bin_exact_failures_curr);
        temperature_bin_shadow_audits_used.push_back(temperature_bin_shadow_audits_curr);
        temperature_bin_shadow_selected_decision_mismatches_used.push_back(temperature_bin_shadow_selected_decision_mismatches_curr);
        temperature_bin_shadow_exact_failures_used.push_back(temperature_bin_shadow_exact_failures_curr);
        temperature_bin_shadow_selected_source_audits_used.push_back(temperature_bin_shadow_selected_source_audits_curr);
        temperature_bin_shadow_selected_source_exact_failures_used.push_back(temperature_bin_shadow_selected_source_exact_failures_curr);

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
                           aux_calibration_points_used,
                           primary_local_calibration_lookups_used,
                           primary_side_calibration_lookups_used,
                           primary_global_calibration_lookups_used,
                           primary_insufficient_calibration_lookups_used,
                           aux_local_calibration_lookups_used,
                           aux_side_calibration_lookups_used,
                           aux_global_calibration_lookups_used,
                           aux_insufficient_calibration_lookups_used,
                           mean_primary_ci_width_used,
                           mean_aux_ci_width_used, mean_selected_ci_width_used,
                           mean_scanned_aux_levels_used, shadow_audits_used,
                           shadow_exact_failures_used, shadow_primary_decision_mismatches_used,
                           shadow_selected_decision_mismatches_used, shadow_selected_interval_hits_used,
                           shadow_full_scan_better_swaps_used, mean_shadow_selected_regret_used,
                           mean_shadow_primary_abs_delta_error_used,
                           mean_shadow_primary_delta_bias_used,
                           mean_shadow_primary_p_bias_used,
                           mean_shadow_selected_abs_delta_error_used,
                           mean_shadow_selected_abs_p_error_used,
                           mean_shadow_selected_delta_bias_used,
                           mean_shadow_selected_p_bias_used,
                           shadow_selected_source_audits_used,
                           shadow_selected_source_exact_failures_used,
                           shadow_selected_source_decision_mismatches_used,
                           shadow_selected_source_interval_hits_used,
                           shadow_selected_source_full_scan_better_swaps_used,
                           mean_shadow_selected_source_regret_used,
                           mean_shadow_selected_source_abs_delta_error_used,
                           mean_shadow_selected_source_abs_p_error_used,
                           mean_shadow_selected_source_delta_bias_used,
                           mean_shadow_selected_source_p_bias_used,
                           temperature_bin_total_swaps_used, temperature_bin_accepted_swaps_used,
                           temperature_bin_primary_resolved_swaps_used, temperature_bin_aux_resolved_swaps_used,
                           temperature_bin_exact_evals_used, temperature_bin_exact_unavailable_swaps_used,
                           temperature_bin_exact_failures_used,
                           temperature_bin_shadow_audits_used,
                           temperature_bin_shadow_selected_decision_mismatches_used,
                           temperature_bin_shadow_exact_failures_used,
                           temperature_bin_shadow_selected_source_audits_used,
                           temperature_bin_shadow_selected_source_exact_failures_used);
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

size_t subsample::aux_level_dim(size_t level) const {
    return test->aux_statistic_level_dim(level);
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

size_t subsample::primary_local_calibration_lookups(size_t k) const {
    return primary_local_calibration_lookups_used.at(k);
}

size_t subsample::primary_side_calibration_lookups(size_t k) const {
    return primary_side_calibration_lookups_used.at(k);
}

size_t subsample::primary_global_calibration_lookups(size_t k) const {
    return primary_global_calibration_lookups_used.at(k);
}

size_t subsample::primary_insufficient_calibration_lookups(size_t k) const {
    return primary_insufficient_calibration_lookups_used.at(k);
}

size_t subsample::aux_local_calibration_lookups(size_t k) const {
    return aux_local_calibration_lookups_used.at(k);
}

size_t subsample::aux_side_calibration_lookups(size_t k) const {
    return aux_side_calibration_lookups_used.at(k);
}

size_t subsample::aux_global_calibration_lookups(size_t k) const {
    return aux_global_calibration_lookups_used.at(k);
}

size_t subsample::aux_insufficient_calibration_lookups(size_t k) const {
    return aux_insufficient_calibration_lookups_used.at(k);
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

double subsample::mean_scanned_aux_levels(size_t k) const {
    return mean_scanned_aux_levels_used.at(k);
}

size_t subsample::shadow_audits(size_t k) const {
    return shadow_audits_used.at(k);
}

size_t subsample::shadow_exact_failures(size_t k) const {
    return shadow_exact_failures_used.at(k);
}

size_t subsample::shadow_primary_decision_mismatches(size_t k) const {
    return shadow_primary_decision_mismatches_used.at(k);
}

size_t subsample::shadow_selected_decision_mismatches(size_t k) const {
    return shadow_selected_decision_mismatches_used.at(k);
}

size_t subsample::shadow_selected_interval_hits(size_t k) const {
    return shadow_selected_interval_hits_used.at(k);
}

size_t subsample::shadow_full_scan_better_swaps(size_t k) const {
    return shadow_full_scan_better_swaps_used.at(k);
}

double subsample::mean_shadow_selected_regret(size_t k) const {
    return mean_shadow_selected_regret_used.at(k);
}

double subsample::mean_shadow_primary_abs_delta_error(size_t k) const {
    return mean_shadow_primary_abs_delta_error_used.at(k);
}

double subsample::mean_shadow_primary_delta_bias(size_t k) const {
    return mean_shadow_primary_delta_bias_used.at(k);
}

double subsample::mean_shadow_primary_p_bias(size_t k) const {
    return mean_shadow_primary_p_bias_used.at(k);
}

double subsample::mean_shadow_selected_abs_delta_error(size_t k) const {
    return mean_shadow_selected_abs_delta_error_used.at(k);
}

double subsample::mean_shadow_selected_abs_p_error(size_t k) const {
    return mean_shadow_selected_abs_p_error_used.at(k);
}

double subsample::mean_shadow_selected_delta_bias(size_t k) const {
    return mean_shadow_selected_delta_bias_used.at(k);
}

double subsample::mean_shadow_selected_p_bias(size_t k) const {
    return mean_shadow_selected_p_bias_used.at(k);
}

size_t subsample::shadow_selected_source_audits(size_t k, size_t source) const {
    return shadow_selected_source_audits_used.at(k).at(source);
}

size_t subsample::shadow_selected_source_exact_failures(size_t k, size_t source) const {
    return shadow_selected_source_exact_failures_used.at(k).at(source);
}

size_t subsample::shadow_selected_source_decision_mismatches(size_t k, size_t source) const {
    return shadow_selected_source_decision_mismatches_used.at(k).at(source);
}

size_t subsample::shadow_selected_source_interval_hits(size_t k, size_t source) const {
    return shadow_selected_source_interval_hits_used.at(k).at(source);
}

size_t subsample::shadow_selected_source_full_scan_better_swaps(size_t k, size_t source) const {
    return shadow_selected_source_full_scan_better_swaps_used.at(k).at(source);
}

double subsample::mean_shadow_selected_source_regret(size_t k, size_t source) const {
    return mean_shadow_selected_source_regret_used.at(k).at(source);
}

double subsample::mean_shadow_selected_source_abs_delta_error(size_t k, size_t source) const {
    return mean_shadow_selected_source_abs_delta_error_used.at(k).at(source);
}

double subsample::mean_shadow_selected_source_abs_p_error(size_t k, size_t source) const {
    return mean_shadow_selected_source_abs_p_error_used.at(k).at(source);
}

double subsample::mean_shadow_selected_source_delta_bias(size_t k, size_t source) const {
    return mean_shadow_selected_source_delta_bias_used.at(k).at(source);
}

double subsample::mean_shadow_selected_source_p_bias(size_t k, size_t source) const {
    return mean_shadow_selected_source_p_bias_used.at(k).at(source);
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

size_t subsample::temperature_bin_shadow_audits(size_t k, size_t bin) const {
    return temperature_bin_shadow_audits_used.at(k).at(bin);
}

size_t subsample::temperature_bin_shadow_selected_decision_mismatches(size_t k, size_t bin) const {
    return temperature_bin_shadow_selected_decision_mismatches_used.at(k).at(bin);
}

size_t subsample::temperature_bin_shadow_exact_failures(size_t k, size_t bin) const {
    return temperature_bin_shadow_exact_failures_used.at(k).at(bin);
}

size_t subsample::temperature_bin_shadow_selected_source_audits(size_t k, size_t bin, size_t source) const {
    const size_t index = bin * CALIBRATION_LOOKUP_SOURCE_COUNT + source;
    return temperature_bin_shadow_selected_source_audits_used.at(k).at(index);
}

size_t subsample::temperature_bin_shadow_selected_source_exact_failures(size_t k, size_t bin, size_t source) const {
    const size_t index = bin * CALIBRATION_LOOKUP_SOURCE_COUNT + source;
    return temperature_bin_shadow_selected_source_exact_failures_used.at(k).at(index);
}

}