#include <Rcpp.h>
#include <vector>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <thread>

#include "include/qchisq.h"
#include "include/matching.h"
#include "include/hw.h"

using namespace Rcpp;
using std::vector;

// [[Rcpp::export]]
LogicalVector quality_control_impl(const IntegerMatrix& case_counts, const NumericVector& maf,
                                   const IntegerVector& mac, const NumericVector& chi2boundary) {
    int n = case_counts.nrow();
    double af = maf[0];
    int ac = mac[0];
    double boundary = chi2boundary[0];
    LogicalVector ret(n);
    for (int i = 0; i < n; i++) {
        int homref = case_counts(i, 0);
        int het = case_counts(i, 1);
        int hom = case_counts(i, 2);
        ret[i] = matching::check_counts(homref, het, hom, af, ac, boundary);
    }
    return ret;
}

namespace {

template<typename F, typename T, typename D>
std::shared_ptr<T> r_to_cpp_impl(const F& matrix, D default_value) {
    std::shared_ptr<T> eigen_matrix = std::make_shared<T>(matrix.nrow(), matrix.ncol());
    for (int i = 0; i < matrix.nrow(); i++) {
        for (int j = 0; j < matrix.ncol(); j++) {
            if (matrix(i, j) == Rcpp::NA) {
                eigen_matrix->operator()(i, j) = default_value;
            } else {
                eigen_matrix->operator()(i, j) = matrix(i, j);
            }
        }
    }
    return eigen_matrix;
}

std::vector<matching::Counts> matrix_to_counts(const Eigen::MatrixXi& matrix) {
    std::vector<matching::Counts> ret;
    if (matrix.cols() != 3) {
        throw std::invalid_argument("Case counts matrix must have three columns");
    }
    for(int i = 0; i < matrix.rows(); i++) {
        auto& row = matrix.row(i);
        ret.emplace_back(row[0], row[1], row[2]);
    }

    return ret;
}

}

std::shared_ptr<Eigen::MatrixXi> r_to_cpp(const IntegerMatrix& matrix) {
    return r_to_cpp_impl<IntegerMatrix, Eigen::MatrixXi, int>(matrix, -1);
}

std::shared_ptr<Eigen::MatrixXd> r_to_cpp(const NumericMatrix& matrix) {
    return r_to_cpp_impl<NumericMatrix, Eigen::MatrixXd, double>(matrix, std::numeric_limits<double>::quiet_NaN());
}

mvn::Vector r_to_cpp(const NumericVector& vector) {
    mvn::Vector eigen_vector(vector.size());
    for (int i = 0; i < vector.size(); i++) {
        eigen_vector(i) = vector(i);
    }
    return eigen_vector;
}


std::vector<std::vector<int>> r_to_cpp_vector(IntegerMatrix& matrix) {
    std::vector<std::vector<int>> ret(matrix.nrow());
    for (int i = 0; i < matrix.nrow(); i++) {
        ret[i].resize(matrix.ncol());
        for (int j = 0; j < matrix.ncol(); j++) {
            ret[i][j] = (matrix(i, j) == Rcpp::NA) ? -1 : matrix(i, j);
        }
    }
    return ret;
}

// [[Rcpp::export]]
List subsample_mvn(NumericMatrix& matrix, IntegerVector size, NumericVector& mean, NumericMatrix& cov) {
    if (matrix.nrow() == 0 || matrix.ncol() == 0) {
        stop("Matrix must be non-empty.");
    }
    if (size.size() != 1) {
        stop("size must contain exactly one value.");
    }
    if (mean.size() != matrix.nrow()) {
        stop("mean must have length equal to nrow(matrix).");
    }
    if (cov.nrow() != matrix.nrow() || cov.ncol() != matrix.nrow()) {
        stop("cov must be a square matrix with nrow(matrix) rows.");
    }
    if (size[0] < matrix.nrow() + 1 || size[0] > matrix.ncol()) {
        stop("Requested subsample size must be between nrow(matrix) + 1 and ncol(matrix).");
    }

    std::vector<int> clusters(matrix.ncol());
    std::iota(clusters.begin(), clusters.end(), 0);
    mvn::Clustering clustering(clusters);
    mvn::subsample annealing(r_to_cpp(matrix), clustering, r_to_cpp(mean), *r_to_cpp(cov), false, 1024);
    annealing.run(1'000'000, 4, 1.0, 0.99995, std::thread::hardware_concurrency(), size[0], size[0], 1);

    List ret;
    auto points = annealing.get_solution(0);
    IntegerVector subset(points.begin(), points.end());
    ret["points"] = subset + 1;

    return ret;
}

// [[Rcpp::export]]
List select_controls_cpp(IntegerMatrix& gmatrix,
                     NumericMatrix& gmatrix_rs,
                     NumericVector& mean, NumericMatrix& directions,
                     IntegerMatrix& cc, IntegerVector& clustering,
                     NumericVector& chi2fn,
                     double min_lambda, double lb_lambda,
                     double max_lambda, double ub_lambda,
                     int min, int max, int step,
                     int sa_iterations, double min_call_rate,
                     std::string method, int n_features,
                     double ci_width_threshold, double calibration_quantile,
                     int min_calibration_samples, int max_calibration_history) {
    vector<double> precomputed_chi(chi2fn.begin(), chi2fn.end());
    qchi2 q(precomputed_chi);

    auto gmatrix_counts = r_to_cpp_vector(gmatrix);
    auto case_counts = r_to_cpp(cc);
    auto principal_directions = r_to_cpp(directions);
    auto gm_rs = r_to_cpp(gmatrix_rs);

    vector<int> clust_vec(clustering.begin(), clustering.end());

    int min_controls = min;
    int max_controls = max;
    int step_clusters = step;
    int iterations = sa_iterations;
    double mcr = min_call_rate;
    int features = n_features;
    int calibration_min = min_calibration_samples;
    int calibration_history = max_calibration_history;
    bool use_nystrom = false;
    bool use_hybrid = false;
    size_t rff_features = std::max(64, features / 4);
    if (method == "exact") {
        use_nystrom = false;
        use_hybrid = false;
    } else if (method == "nystrom") {
        use_nystrom = true;
        use_hybrid = false;
    } else if (method == "hybrid") {
        use_nystrom = true;
        use_hybrid = true;
        rff_features = std::min<size_t>(20000, std::max<size_t>(1024, 16 * (size_t)features));
    } else {
        stop("Unsupported method. Use 'exact', 'nystrom', or 'hybrid'.");
    }
    mvn::Clustering cl(clust_vec);

    matching::matching matcher(std::move(gmatrix_counts), gm_rs, cl);
    matcher.set_qchi_sq_function(q.function());
    matcher.set_soft_threshold({lb_lambda, ub_lambda});
    matcher.set_hard_threshold({min_lambda, max_lambda});
    matcher.process_mvn(*principal_directions, r_to_cpp(mean), std::thread::hardware_concurrency(),
                        min_controls, max_controls, step_clusters, iterations, use_nystrom, (size_t)features,
                        use_hybrid, rff_features, 42u,
                        ci_width_threshold, calibration_quantile,
                        (size_t)calibration_min, (size_t)calibration_history);
    matcher.set_interrupts_checker([]() { Rcpp::checkUserInterrupt(); });

    auto result = matcher.match(matrix_to_counts(*case_counts), min_controls, mcr);

    List ret;
    NumericVector lambda(result.lambdas.begin(), result.lambdas.end());
    NumericVector pvals(result.pvals.begin(), result.pvals.end());
    NumericVector stats(result.statistics.begin(), result.statistics.end());
    IntegerVector names(result.lambda_i.begin(), result.lambda_i.end());
    IntegerVector pvals_num(result.pvals_num.begin(), result.pvals_num.end());
    IntegerVector optimal_controls(result.optimal_prefix.begin(), result.optimal_prefix.end());
    NumericVector optimal_lambda = {result.optimal_lambda};

    const size_t n_sa = matcher.sa_solutions();
    IntegerVector sa_sizes((R_xlen_t)n_sa);
    IntegerVector sa_total_swaps((R_xlen_t)n_sa);
    IntegerVector sa_uncertain_swaps((R_xlen_t)n_sa);
    IntegerVector sa_ci_resolved_swaps((R_xlen_t)n_sa);
    IntegerVector sa_primary_resolved_swaps((R_xlen_t)n_sa);
    IntegerVector sa_exact_unavailable_swaps((R_xlen_t)n_sa);
    IntegerVector sa_exact_evals((R_xlen_t)n_sa);
    IntegerVector sa_exact_evals_on_improving((R_xlen_t)n_sa);
    IntegerVector sa_exact_evals_on_worsening((R_xlen_t)n_sa);
    IntegerVector sa_exact_failures((R_xlen_t)n_sa);
    IntegerVector sa_primary_calibration_points((R_xlen_t)n_sa);
    IntegerVector sa_aux_calibration_points((R_xlen_t)n_sa);
    IntegerVector sa_primary_local_calibration_lookups((R_xlen_t)n_sa);
    IntegerVector sa_primary_side_calibration_lookups((R_xlen_t)n_sa);
    IntegerVector sa_primary_global_calibration_lookups((R_xlen_t)n_sa);
    IntegerVector sa_primary_insufficient_calibration_lookups((R_xlen_t)n_sa);
    IntegerVector sa_aux_local_calibration_lookups((R_xlen_t)n_sa);
    IntegerVector sa_aux_side_calibration_lookups((R_xlen_t)n_sa);
    IntegerVector sa_aux_global_calibration_lookups((R_xlen_t)n_sa);
    IntegerVector sa_aux_insufficient_calibration_lookups((R_xlen_t)n_sa);
    NumericVector sa_mean_primary_ci_width((R_xlen_t)n_sa);
    NumericVector sa_mean_aux_ci_width((R_xlen_t)n_sa);
    NumericVector sa_mean_selected_ci_width((R_xlen_t)n_sa);
    NumericVector sa_mean_scanned_aux_levels((R_xlen_t)n_sa);
    IntegerVector sa_shadow_audits((R_xlen_t)n_sa);
    IntegerVector sa_shadow_exact_failures((R_xlen_t)n_sa);
    IntegerVector sa_shadow_primary_decision_mismatches((R_xlen_t)n_sa);
    IntegerVector sa_shadow_selected_decision_mismatches((R_xlen_t)n_sa);
    IntegerVector sa_shadow_selected_interval_hits((R_xlen_t)n_sa);
    NumericVector sa_shadow_selected_interval_coverage((R_xlen_t)n_sa);
    IntegerVector sa_shadow_full_scan_better_swaps((R_xlen_t)n_sa);
    NumericVector sa_mean_shadow_selected_regret((R_xlen_t)n_sa);
    NumericVector sa_mean_shadow_primary_abs_delta_error((R_xlen_t)n_sa);
    NumericVector sa_mean_shadow_primary_delta_bias((R_xlen_t)n_sa);
    NumericVector sa_mean_shadow_primary_p_bias((R_xlen_t)n_sa);
    NumericVector sa_mean_shadow_selected_abs_delta_error((R_xlen_t)n_sa);
    NumericVector sa_mean_shadow_selected_abs_p_error((R_xlen_t)n_sa);
    NumericVector sa_mean_shadow_selected_delta_bias((R_xlen_t)n_sa);
    NumericVector sa_mean_shadow_selected_p_bias((R_xlen_t)n_sa);
    CharacterVector sa_names((R_xlen_t)n_sa);
    const size_t n_aux_levels = matcher.sa_aux_ladder_levels();
    List sa_aux_level_resolved((R_xlen_t)n_aux_levels);
    CharacterVector sa_aux_level_names((R_xlen_t)n_aux_levels);
    for (size_t level = 0; level < n_aux_levels; ++level) {
        IntegerVector level_counts((R_xlen_t)n_sa);
        for (size_t i = 0; i < n_sa; ++i) {
            level_counts[(R_xlen_t)i] = (int)matcher.sa_aux_level_resolved_swaps(i, level);
        }
        sa_aux_level_resolved[(R_xlen_t)level] = level_counts;
        sa_aux_level_names[(R_xlen_t)level] = "f" + std::to_string(matcher.sa_aux_level_dim(level));
    }
    sa_aux_level_resolved.attr("names") = sa_aux_level_names;

    const size_t n_shadow_sources = 4;
    const char* shadow_source_names_raw[] = {"local", "side", "global", "insufficient"};
    CharacterVector sa_shadow_source_names((R_xlen_t)n_shadow_sources);
    for (size_t source = 0; source < n_shadow_sources; ++source) {
        sa_shadow_source_names[(R_xlen_t)source] = shadow_source_names_raw[source];
    }

    const size_t n_temp_bins = matcher.sa_temperature_bins();
    List sa_temperature_bins((R_xlen_t)n_temp_bins);
    CharacterVector sa_temperature_bin_names((R_xlen_t)n_temp_bins);
    for (size_t i = 0; i < n_sa; ++i) {
        const int size = (int)matcher.sa_solution_size(i);
        sa_sizes[(R_xlen_t)i] = size;
        sa_total_swaps[(R_xlen_t)i] = (int)matcher.sa_total_swaps(i);
        sa_uncertain_swaps[(R_xlen_t)i] = (int)matcher.sa_uncertain_swaps(i);
        sa_ci_resolved_swaps[(R_xlen_t)i] = (int)matcher.sa_ci_resolved_swaps(i);
        sa_primary_resolved_swaps[(R_xlen_t)i] = (int)matcher.sa_primary_resolved_swaps(i);
        sa_exact_unavailable_swaps[(R_xlen_t)i] = (int)matcher.sa_exact_unavailable_swaps(i);
        sa_exact_evals[(R_xlen_t)i] = (int)matcher.sa_exact_evals(i);
        sa_exact_evals_on_improving[(R_xlen_t)i] = (int)matcher.sa_exact_evals_on_improving(i);
        sa_exact_evals_on_worsening[(R_xlen_t)i] = (int)matcher.sa_exact_evals_on_worsening(i);
        sa_exact_failures[(R_xlen_t)i] = (int)matcher.sa_exact_failures(i);
        sa_primary_calibration_points[(R_xlen_t)i] = (int)matcher.sa_primary_calibration_points(i);
        sa_aux_calibration_points[(R_xlen_t)i] = (int)matcher.sa_aux_calibration_points(i);
        sa_primary_local_calibration_lookups[(R_xlen_t)i] = (int)matcher.sa_primary_local_calibration_lookups(i);
        sa_primary_side_calibration_lookups[(R_xlen_t)i] = (int)matcher.sa_primary_side_calibration_lookups(i);
        sa_primary_global_calibration_lookups[(R_xlen_t)i] = (int)matcher.sa_primary_global_calibration_lookups(i);
        sa_primary_insufficient_calibration_lookups[(R_xlen_t)i] = (int)matcher.sa_primary_insufficient_calibration_lookups(i);
        sa_aux_local_calibration_lookups[(R_xlen_t)i] = (int)matcher.sa_aux_local_calibration_lookups(i);
        sa_aux_side_calibration_lookups[(R_xlen_t)i] = (int)matcher.sa_aux_side_calibration_lookups(i);
        sa_aux_global_calibration_lookups[(R_xlen_t)i] = (int)matcher.sa_aux_global_calibration_lookups(i);
        sa_aux_insufficient_calibration_lookups[(R_xlen_t)i] = (int)matcher.sa_aux_insufficient_calibration_lookups(i);
        sa_mean_primary_ci_width[(R_xlen_t)i] = matcher.sa_mean_primary_ci_width(i);
        sa_mean_aux_ci_width[(R_xlen_t)i] = matcher.sa_mean_aux_ci_width(i);
        sa_mean_selected_ci_width[(R_xlen_t)i] = matcher.sa_mean_selected_ci_width(i);
        sa_mean_scanned_aux_levels[(R_xlen_t)i] = matcher.sa_mean_scanned_aux_levels(i);
        sa_shadow_audits[(R_xlen_t)i] = (int)matcher.sa_shadow_audits(i);
        sa_shadow_exact_failures[(R_xlen_t)i] = (int)matcher.sa_shadow_exact_failures(i);
        sa_shadow_primary_decision_mismatches[(R_xlen_t)i] = (int)matcher.sa_shadow_primary_decision_mismatches(i);
        sa_shadow_selected_decision_mismatches[(R_xlen_t)i] = (int)matcher.sa_shadow_selected_decision_mismatches(i);
        sa_shadow_selected_interval_hits[(R_xlen_t)i] = (int)matcher.sa_shadow_selected_interval_hits(i);
        const size_t shadow_successes = matcher.sa_shadow_audits(i) - matcher.sa_shadow_exact_failures(i);
        sa_shadow_selected_interval_coverage[(R_xlen_t)i] = shadow_successes == 0
            ? NA_REAL
            : static_cast<double>(matcher.sa_shadow_selected_interval_hits(i)) /
                  static_cast<double>(shadow_successes);
        sa_shadow_full_scan_better_swaps[(R_xlen_t)i] = (int)matcher.sa_shadow_full_scan_better_swaps(i);
        sa_mean_shadow_selected_regret[(R_xlen_t)i] = matcher.sa_mean_shadow_selected_regret(i);
        sa_mean_shadow_primary_abs_delta_error[(R_xlen_t)i] = matcher.sa_mean_shadow_primary_abs_delta_error(i);
        sa_mean_shadow_primary_delta_bias[(R_xlen_t)i] = matcher.sa_mean_shadow_primary_delta_bias(i);
        sa_mean_shadow_primary_p_bias[(R_xlen_t)i] = matcher.sa_mean_shadow_primary_p_bias(i);
        sa_mean_shadow_selected_abs_delta_error[(R_xlen_t)i] = matcher.sa_mean_shadow_selected_abs_delta_error(i);
        sa_mean_shadow_selected_abs_p_error[(R_xlen_t)i] = matcher.sa_mean_shadow_selected_abs_p_error(i);
        sa_mean_shadow_selected_delta_bias[(R_xlen_t)i] = matcher.sa_mean_shadow_selected_delta_bias(i);
        sa_mean_shadow_selected_p_bias[(R_xlen_t)i] = matcher.sa_mean_shadow_selected_p_bias(i);
        sa_names[(R_xlen_t)i] = std::to_string(size);
    }
    for (size_t level = 0; level < n_aux_levels; ++level) {
        IntegerVector level_counts = sa_aux_level_resolved[(R_xlen_t)level];
        level_counts.attr("names") = sa_names;
        sa_aux_level_resolved[(R_xlen_t)level] = level_counts;
    }
    for (size_t bin = 0; bin < n_temp_bins; ++bin) {
        List bin_stats = List::create(
            Named("total_swaps") = IntegerVector((R_xlen_t)n_sa),
            Named("accepted_swaps") = IntegerVector((R_xlen_t)n_sa),
            Named("primary_resolved_swaps") = IntegerVector((R_xlen_t)n_sa),
            Named("aux_resolved_swaps") = IntegerVector((R_xlen_t)n_sa),
            Named("exact_evals") = IntegerVector((R_xlen_t)n_sa),
            Named("exact_unavailable_swaps") = IntegerVector((R_xlen_t)n_sa),
            Named("exact_failures") = IntegerVector((R_xlen_t)n_sa),
            Named("shadow_audits") = IntegerVector((R_xlen_t)n_sa),
            Named("shadow_selected_decision_mismatches") = IntegerVector((R_xlen_t)n_sa),
            Named("shadow_exact_failures") = IntegerVector((R_xlen_t)n_sa),
            Named("shadow_selected_by_source") = List::create()
        );
        IntegerVector total_swaps = bin_stats["total_swaps"];
        IntegerVector accepted_swaps = bin_stats["accepted_swaps"];
        IntegerVector primary_resolved_swaps = bin_stats["primary_resolved_swaps"];
        IntegerVector aux_resolved_swaps = bin_stats["aux_resolved_swaps"];
        IntegerVector exact_evals = bin_stats["exact_evals"];
        IntegerVector exact_unavailable_swaps = bin_stats["exact_unavailable_swaps"];
        IntegerVector exact_failures = bin_stats["exact_failures"];
        IntegerVector shadow_audits = bin_stats["shadow_audits"];
        IntegerVector shadow_selected_decision_mismatches = bin_stats["shadow_selected_decision_mismatches"];
        IntegerVector shadow_exact_failures = bin_stats["shadow_exact_failures"];
        for (size_t i = 0; i < n_sa; ++i) {
            total_swaps[(R_xlen_t)i] = (int)matcher.sa_temperature_bin_total_swaps(i, bin);
            accepted_swaps[(R_xlen_t)i] = (int)matcher.sa_temperature_bin_accepted_swaps(i, bin);
            primary_resolved_swaps[(R_xlen_t)i] = (int)matcher.sa_temperature_bin_primary_resolved_swaps(i, bin);
            aux_resolved_swaps[(R_xlen_t)i] = (int)matcher.sa_temperature_bin_aux_resolved_swaps(i, bin);
            exact_evals[(R_xlen_t)i] = (int)matcher.sa_temperature_bin_exact_evals(i, bin);
            exact_unavailable_swaps[(R_xlen_t)i] = (int)matcher.sa_temperature_bin_exact_unavailable_swaps(i, bin);
            exact_failures[(R_xlen_t)i] = (int)matcher.sa_temperature_bin_exact_failures(i, bin);
            shadow_audits[(R_xlen_t)i] = (int)matcher.sa_temperature_bin_shadow_audits(i, bin);
            shadow_selected_decision_mismatches[(R_xlen_t)i] =
                (int)matcher.sa_temperature_bin_shadow_selected_decision_mismatches(i, bin);
            shadow_exact_failures[(R_xlen_t)i] = (int)matcher.sa_temperature_bin_shadow_exact_failures(i, bin);
        }
        total_swaps.attr("names") = sa_names;
        accepted_swaps.attr("names") = sa_names;
        primary_resolved_swaps.attr("names") = sa_names;
        aux_resolved_swaps.attr("names") = sa_names;
        exact_evals.attr("names") = sa_names;
        exact_unavailable_swaps.attr("names") = sa_names;
        exact_failures.attr("names") = sa_names;
        shadow_audits.attr("names") = sa_names;
        shadow_selected_decision_mismatches.attr("names") = sa_names;
        shadow_exact_failures.attr("names") = sa_names;

        List shadow_selected_by_source((R_xlen_t)n_shadow_sources);
        for (size_t source = 0; source < n_shadow_sources; ++source) {
            IntegerVector source_audits((R_xlen_t)n_sa);
            IntegerVector source_exact_failures((R_xlen_t)n_sa);
            IntegerVector source_successful_audits((R_xlen_t)n_sa);
            for (size_t i = 0; i < n_sa; ++i) {
                const size_t audits = matcher.sa_temperature_bin_shadow_selected_source_audits(i, bin, source);
                const size_t exact_failures_by_source = matcher.sa_temperature_bin_shadow_selected_source_exact_failures(i, bin, source);
                source_audits[(R_xlen_t)i] = (int)audits;
                source_exact_failures[(R_xlen_t)i] = (int)exact_failures_by_source;
                source_successful_audits[(R_xlen_t)i] = (int)(audits - exact_failures_by_source);
            }
            source_audits.attr("names") = sa_names;
            source_exact_failures.attr("names") = sa_names;
            source_successful_audits.attr("names") = sa_names;
            shadow_selected_by_source[(R_xlen_t)source] = List::create(
                Named("audits") = source_audits,
                Named("exact_failures") = source_exact_failures,
                Named("successful_audits") = source_successful_audits
            );
        }
        shadow_selected_by_source.attr("names") = sa_shadow_source_names;
        bin_stats["total_swaps"] = total_swaps;
        bin_stats["accepted_swaps"] = accepted_swaps;
        bin_stats["primary_resolved_swaps"] = primary_resolved_swaps;
        bin_stats["aux_resolved_swaps"] = aux_resolved_swaps;
        bin_stats["exact_evals"] = exact_evals;
        bin_stats["exact_unavailable_swaps"] = exact_unavailable_swaps;
        bin_stats["exact_failures"] = exact_failures;
        bin_stats["shadow_audits"] = shadow_audits;
        bin_stats["shadow_selected_decision_mismatches"] = shadow_selected_decision_mismatches;
        bin_stats["shadow_exact_failures"] = shadow_exact_failures;
        bin_stats["shadow_selected_by_source"] = shadow_selected_by_source;
        sa_temperature_bins[(R_xlen_t)bin] = bin_stats;
        sa_temperature_bin_names[(R_xlen_t)bin] = "bin_" + std::to_string(bin + 1);
    }
    sa_temperature_bins.attr("names") = sa_temperature_bin_names;

    List sa_shadow_selected_by_source((R_xlen_t)n_shadow_sources);
    for (size_t source = 0; source < n_shadow_sources; ++source) {
        IntegerVector source_audits((R_xlen_t)n_sa);
        IntegerVector source_exact_failures((R_xlen_t)n_sa);
        IntegerVector source_decision_mismatches((R_xlen_t)n_sa);
        IntegerVector source_interval_hits((R_xlen_t)n_sa);
        IntegerVector source_full_scan_better_swaps((R_xlen_t)n_sa);
        NumericVector source_interval_coverage((R_xlen_t)n_sa);
        NumericVector source_decision_mismatch_rate((R_xlen_t)n_sa);
        NumericVector source_full_scan_better_rate((R_xlen_t)n_sa);
        NumericVector source_mean_regret((R_xlen_t)n_sa);
        NumericVector source_mean_abs_delta_error((R_xlen_t)n_sa);
        NumericVector source_mean_abs_p_error((R_xlen_t)n_sa);
        NumericVector source_mean_delta_bias((R_xlen_t)n_sa);
        NumericVector source_mean_p_bias((R_xlen_t)n_sa);
        for (size_t i = 0; i < n_sa; ++i) {
            const size_t audits = matcher.sa_shadow_selected_source_audits(i, source);
            const size_t exact_failures = matcher.sa_shadow_selected_source_exact_failures(i, source);
            const size_t successes = audits - exact_failures;
            const size_t decision_mismatches = matcher.sa_shadow_selected_source_decision_mismatches(i, source);
            const size_t interval_hits = matcher.sa_shadow_selected_source_interval_hits(i, source);
            const size_t full_scan_better_swaps = matcher.sa_shadow_selected_source_full_scan_better_swaps(i, source);

            source_audits[(R_xlen_t)i] = (int)audits;
            source_exact_failures[(R_xlen_t)i] = (int)exact_failures;
            source_decision_mismatches[(R_xlen_t)i] = (int)decision_mismatches;
            source_interval_hits[(R_xlen_t)i] = (int)interval_hits;
            source_full_scan_better_swaps[(R_xlen_t)i] = (int)full_scan_better_swaps;
            source_interval_coverage[(R_xlen_t)i] = successes == 0
                ? NA_REAL
                : static_cast<double>(interval_hits) / static_cast<double>(successes);
            source_decision_mismatch_rate[(R_xlen_t)i] = successes == 0
                ? NA_REAL
                : static_cast<double>(decision_mismatches) / static_cast<double>(successes);
            source_full_scan_better_rate[(R_xlen_t)i] = successes == 0
                ? NA_REAL
                : static_cast<double>(full_scan_better_swaps) / static_cast<double>(successes);
            source_mean_regret[(R_xlen_t)i] = matcher.sa_mean_shadow_selected_source_regret(i, source);
            source_mean_abs_delta_error[(R_xlen_t)i] = matcher.sa_mean_shadow_selected_source_abs_delta_error(i, source);
            source_mean_abs_p_error[(R_xlen_t)i] = matcher.sa_mean_shadow_selected_source_abs_p_error(i, source);
            source_mean_delta_bias[(R_xlen_t)i] = matcher.sa_mean_shadow_selected_source_delta_bias(i, source);
            source_mean_p_bias[(R_xlen_t)i] = matcher.sa_mean_shadow_selected_source_p_bias(i, source);
        }

        source_audits.attr("names") = sa_names;
        source_exact_failures.attr("names") = sa_names;
        source_decision_mismatches.attr("names") = sa_names;
        source_interval_hits.attr("names") = sa_names;
        source_full_scan_better_swaps.attr("names") = sa_names;
        source_interval_coverage.attr("names") = sa_names;
        source_decision_mismatch_rate.attr("names") = sa_names;
        source_full_scan_better_rate.attr("names") = sa_names;
        source_mean_regret.attr("names") = sa_names;
        source_mean_abs_delta_error.attr("names") = sa_names;
        source_mean_abs_p_error.attr("names") = sa_names;
        source_mean_delta_bias.attr("names") = sa_names;
        source_mean_p_bias.attr("names") = sa_names;

        sa_shadow_selected_by_source[(R_xlen_t)source] = List::create(
            Named("audits") = source_audits,
            Named("exact_failures") = source_exact_failures,
            Named("decision_mismatches") = source_decision_mismatches,
            Named("interval_hits") = source_interval_hits,
            Named("full_scan_better_swaps") = source_full_scan_better_swaps,
            Named("interval_coverage") = source_interval_coverage,
            Named("decision_mismatch_rate") = source_decision_mismatch_rate,
            Named("full_scan_better_rate") = source_full_scan_better_rate,
            Named("mean_regret") = source_mean_regret,
            Named("mean_abs_delta_error") = source_mean_abs_delta_error,
            Named("mean_abs_p_error") = source_mean_abs_p_error,
            Named("mean_delta_bias") = source_mean_delta_bias,
            Named("mean_p_bias") = source_mean_p_bias
        );
    }
    sa_shadow_selected_by_source.attr("names") = sa_shadow_source_names;

    lambda.attr("names") = names;
    stats.attr("names") = names;
    pvals_num.attr("names") = names;
    sa_total_swaps.attr("names") = sa_names;
    sa_uncertain_swaps.attr("names") = sa_names;
    sa_ci_resolved_swaps.attr("names") = sa_names;
    sa_primary_resolved_swaps.attr("names") = sa_names;
    sa_exact_unavailable_swaps.attr("names") = sa_names;
    sa_exact_evals.attr("names") = sa_names;
    sa_exact_evals_on_improving.attr("names") = sa_names;
    sa_exact_evals_on_worsening.attr("names") = sa_names;
    sa_exact_failures.attr("names") = sa_names;
    sa_primary_calibration_points.attr("names") = sa_names;
    sa_aux_calibration_points.attr("names") = sa_names;
    sa_primary_local_calibration_lookups.attr("names") = sa_names;
    sa_primary_side_calibration_lookups.attr("names") = sa_names;
    sa_primary_global_calibration_lookups.attr("names") = sa_names;
    sa_primary_insufficient_calibration_lookups.attr("names") = sa_names;
    sa_aux_local_calibration_lookups.attr("names") = sa_names;
    sa_aux_side_calibration_lookups.attr("names") = sa_names;
    sa_aux_global_calibration_lookups.attr("names") = sa_names;
    sa_aux_insufficient_calibration_lookups.attr("names") = sa_names;
    sa_mean_primary_ci_width.attr("names") = sa_names;
    sa_mean_aux_ci_width.attr("names") = sa_names;
    sa_mean_selected_ci_width.attr("names") = sa_names;
    sa_mean_scanned_aux_levels.attr("names") = sa_names;
    sa_shadow_audits.attr("names") = sa_names;
    sa_shadow_exact_failures.attr("names") = sa_names;
    sa_shadow_primary_decision_mismatches.attr("names") = sa_names;
    sa_shadow_selected_decision_mismatches.attr("names") = sa_names;
    sa_shadow_selected_interval_hits.attr("names") = sa_names;
    sa_shadow_selected_interval_coverage.attr("names") = sa_names;
    sa_shadow_full_scan_better_swaps.attr("names") = sa_names;
    sa_mean_shadow_selected_regret.attr("names") = sa_names;
    sa_mean_shadow_primary_abs_delta_error.attr("names") = sa_names;
    sa_mean_shadow_primary_delta_bias.attr("names") = sa_names;
    sa_mean_shadow_primary_p_bias.attr("names") = sa_names;
    sa_mean_shadow_selected_abs_delta_error.attr("names") = sa_names;
    sa_mean_shadow_selected_abs_p_error.attr("names") = sa_names;
    sa_mean_shadow_selected_delta_bias.attr("names") = sa_names;
    sa_mean_shadow_selected_p_bias.attr("names") = sa_names;
    ret["lambda"] = lambda;
    ret["optimal_lambda"] = optimal_lambda;
    ret["statistics"] = stats;
    ret["controls"] = optimal_controls + 1;
    ret["pvals"] = pvals;
    ret["snps"] = pvals_num;
    ret["sa_diagnostics"] = List::create(
        Named("subset_sizes") = sa_sizes,
        Named("total_swaps") = sa_total_swaps,
        Named("uncertain_swaps") = sa_uncertain_swaps,
        Named("ci_resolved_swaps") = sa_ci_resolved_swaps,
        Named("primary_resolved_swaps") = sa_primary_resolved_swaps,
        Named("aux_ladder_resolved_swaps") = sa_aux_level_resolved,
        Named("exact_unavailable_swaps") = sa_exact_unavailable_swaps,
        Named("exact_evals") = sa_exact_evals,
        Named("exact_evals_on_improving_swaps") = sa_exact_evals_on_improving,
        Named("exact_evals_on_worsening_swaps") = sa_exact_evals_on_worsening,
        Named("exact_failures") = sa_exact_failures,
        Named("primary_calibration_points") = sa_primary_calibration_points,
        Named("aux_calibration_points") = sa_aux_calibration_points,
        Named("primary_local_calibration_lookups") = sa_primary_local_calibration_lookups,
        Named("primary_side_calibration_lookups") = sa_primary_side_calibration_lookups,
        Named("primary_global_calibration_lookups") = sa_primary_global_calibration_lookups,
        Named("primary_insufficient_calibration_lookups") = sa_primary_insufficient_calibration_lookups,
        Named("aux_local_calibration_lookups") = sa_aux_local_calibration_lookups,
        Named("aux_side_calibration_lookups") = sa_aux_side_calibration_lookups,
        Named("aux_global_calibration_lookups") = sa_aux_global_calibration_lookups,
        Named("aux_insufficient_calibration_lookups") = sa_aux_insufficient_calibration_lookups,
        Named("mean_primary_ci_width") = sa_mean_primary_ci_width,
        Named("mean_aux_ci_width") = sa_mean_aux_ci_width,
        Named("mean_selected_ci_width") = sa_mean_selected_ci_width,
        Named("mean_scanned_aux_levels") = sa_mean_scanned_aux_levels,
        Named("shadow_audits") = sa_shadow_audits,
        Named("shadow_exact_failures") = sa_shadow_exact_failures,
        Named("shadow_primary_decision_mismatches") = sa_shadow_primary_decision_mismatches,
        Named("shadow_selected_decision_mismatches") = sa_shadow_selected_decision_mismatches,
        Named("shadow_selected_interval_hits") = sa_shadow_selected_interval_hits,
        Named("shadow_selected_interval_coverage") = sa_shadow_selected_interval_coverage,
        Named("shadow_full_scan_better_swaps") = sa_shadow_full_scan_better_swaps,
        Named("mean_shadow_selected_regret") = sa_mean_shadow_selected_regret,
        Named("mean_shadow_primary_abs_delta_error") = sa_mean_shadow_primary_abs_delta_error,
        Named("mean_shadow_primary_delta_bias") = sa_mean_shadow_primary_delta_bias,
        Named("mean_shadow_primary_p_bias") = sa_mean_shadow_primary_p_bias,
        Named("mean_shadow_selected_abs_delta_error") = sa_mean_shadow_selected_abs_delta_error,
        Named("mean_shadow_selected_abs_p_error") = sa_mean_shadow_selected_abs_p_error,
        Named("mean_shadow_selected_delta_bias") = sa_mean_shadow_selected_delta_bias,
        Named("mean_shadow_selected_p_bias") = sa_mean_shadow_selected_p_bias,
        Named("shadow_selected_by_source") = sa_shadow_selected_by_source,
        Named("temperature_bins") = sa_temperature_bins
    );
    return ret;
}
