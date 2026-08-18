#include <Rcpp.h>
#include <vector>
#include <RcppEigen.h>
#include <limits>

// [[Rcpp::depends(RcppEigen)]]

#include <thread>

#include "include/qchisq.h"
#include "include/matching.h"
#include "include/subsample.h"
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

constexpr size_t max_stored_cluster_size = std::numeric_limits<std::uint8_t>::max();

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

void validate_cluster_sizes(const mvn::Clustering& clustering) {
    for (size_t cluster = 0; cluster < clustering.size(); cluster++) {
        if (clustering.cluster_size(cluster) > max_stored_cluster_size) {
            Rcpp::stop("Each control cluster must contain at most 255 samples");
        }
    }
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


// `variant_rows` maps every case variant (row of the case count matrix) to the
// row of `matrix` that holds the matching control genotypes, or -1 when the
// case variant is absent from the controls. This lets the caller pass the full
// control genotype matrix without slicing it down to the case variants.
std::vector<std::vector<matching::ClusterCounts>> build_cluster_counts(const IntegerMatrix& matrix,
                                                                       const std::vector<int>& clustering,
                                                                       size_t n_clusters,
                                                                       const std::vector<int>& variant_rows) {
    int n_case_variants = static_cast<int>(variant_rows.size());
    int n_control_variants = matrix.nrow();
    int n_samples = matrix.ncol();
    std::vector<std::vector<matching::ClusterCounts>> counts(n_clusters,
                                                             std::vector<matching::ClusterCounts>(n_case_variants));
    for (int j = 0; j < n_samples; j++) {
        int cluster = clustering[j];
        for (int cv = 0; cv < n_case_variants; cv++) {
            int i = variant_rows[cv];
            if (i < 0 || i >= n_control_variants) {
                continue;
            }
            int value = matrix(i, j);
            if (!(value == Rcpp::NA)) {
                counts[cluster][cv][value] += 1;
            }
        }
    }
    return counts;
}

// [[Rcpp::export]]
List subsample_mvn(NumericMatrix& matrix, IntegerVector size, NumericVector& mean, NumericMatrix& cov,
                   int seed) {
    std::vector<int> clusters(matrix.ncol());
    std::iota(clusters.begin(), clusters.end(), 0);
    mvn::Clustering clustering(clusters);
    mvn::subsample annealing(r_to_cpp(matrix), clustering, r_to_cpp(mean), *r_to_cpp(cov),
                             mvn::PrecomputeConfig{}, static_cast<std::mt19937::result_type>(seed));
    annealing.run(1'000'000, 4, 1.0, 0.99995, std::thread::hardware_concurrency(), size[0], size[0], 1);

    List ret;
    auto points = annealing.get_solution(0);
    IntegerVector subset(points.begin(), points.end());
    ret["points"] = subset + 1;

    return ret;
}

// Second stage of the control selection: the multivariate normality search on its
// own. It knows nothing about genotypes -- it only sees points in some reduced
// space, the clusters they are swapped in as units, and the target normal
// distribution -- and it returns one candidate subset of clusters per target size
// rather than a single answer. Choosing between those candidates is the caller's
// job (see match_controls_cpp).
// [[Rcpp::export]]
List mvn_subsample_clusters_cpp(NumericMatrix& points,
                                IntegerVector& clustering,
                                NumericVector& mean,
                                NumericMatrix& cov,
                                int min, int max, int step,
                                int iterations, int seed,
                                int sa_threads = 0,
                                int exact_precompute_threads = 0,
                                int exact_cluster_tile_size = 32) {
    // The temperature decays geometrically from t0 = 1 down to this floor over the
    // whole iteration budget.
    const double TEMPERATURE_FLOOR = 1e-18;

    vector<int> clust_vec(clustering.begin(), clustering.end());
    mvn::Clustering cl(clust_vec);

    mvn::PrecomputeConfig config;
    if (exact_precompute_threads > 0) {
        config.threads = static_cast<size_t>(exact_precompute_threads);
    }
    if (exact_cluster_tile_size > 0) {
        config.cluster_tile_size = static_cast<size_t>(exact_cluster_tile_size);
    }

    int sa_pool_size = sa_threads > 0 ? sa_threads : static_cast<int>(std::thread::hardware_concurrency());
    if (sa_pool_size <= 0) {
        sa_pool_size = 1;
    }

    auto space = r_to_cpp(points);
    Rcpp::Rcerr << "Starting processing controls space." << std::endl;
    Rcpp::Rcerr << "The size of controls space is " << space->rows() << " by " << space->cols() << std::endl;

    mvn::subsample annealing(space, cl, r_to_cpp(mean), *r_to_cpp(cov), config,
                             static_cast<std::mt19937::result_type>(seed));
    space.reset();
    Rcpp::Rcerr << "Mahalanobis distances have been successfully calculated." << std::endl;

    double c = std::pow(TEMPERATURE_FLOOR, 1.0 / (double)iterations);
    annealing.run(iterations, 4, 1.0, c, sa_pool_size, min, max, step);

    int n_solutions = static_cast<int>(annealing.solutions());
    List clusters(n_solutions);
    NumericVector statistics(n_solutions);
    for (int k = 0; k < n_solutions; k++) {
        auto solution = annealing.get_solution(k);
        IntegerVector ids(solution.begin(), solution.end());
        clusters[k] = ids + 1;
        statistics[k] = annealing.statistic(k);
    }

    List ret;
    ret["clusters"] = clusters;
    ret["statistics"] = statistics;
    return ret;
}

// Third stage: pick the candidate subset whose lambda_GC matches the cases best.
// `candidates` holds one-based cluster ids as produced by mvn_subsample_clusters_cpp,
// and `variant_rows` maps case variants onto rows of `gmatrix` (see build_cluster_counts).
// [[Rcpp::export]]
List match_controls_cpp(IntegerMatrix& gmatrix,
                        IntegerMatrix& cc,
                        IntegerVector& variant_rows,
                        IntegerVector& clustering,
                        List& candidates,
                        NumericVector& statistics,
                        NumericVector& chi2fn,
                        double min_lambda, double lb_lambda,
                        double max_lambda, double ub_lambda,
                        int min_controls, double min_call_rate) {
    vector<double> precomputed_chi(chi2fn.begin(), chi2fn.end());
    qchi2 q(precomputed_chi);

    auto case_counts = r_to_cpp(cc);

    vector<int> clust_vec(clustering.begin(), clustering.end());
    vector<int> variant_rows_vec(variant_rows.begin(), variant_rows.end());

    mvn::Clustering cl(clust_vec);
    validate_cluster_sizes(cl);

    if (candidates.size() != statistics.size()) {
        Rcpp::stop("Every candidate subset must come with its normality statistic");
    }

    std::vector<std::vector<size_t>> candidate_subsets;
    candidate_subsets.reserve(candidates.size());
    for (R_xlen_t k = 0; k < candidates.size(); k++) {
        IntegerVector ids = candidates[k];
        std::vector<size_t> subset;
        subset.reserve(ids.size());
        for (int id: ids) {
            if (id < 1 || static_cast<size_t>(id) > cl.size()) {
                Rcpp::stop("Candidate subsets must reference existing clusters");
            }
            subset.push_back(static_cast<size_t>(id - 1));
        }
        candidate_subsets.push_back(std::move(subset));
    }

    auto cluster_counts = build_cluster_counts(gmatrix, clust_vec, cl.size(), variant_rows_vec);
    matching::matching matcher(std::move(cluster_counts), cl);
    matcher.set_qchi_sq_function(q.function());
    matcher.set_soft_threshold({lb_lambda, ub_lambda});
    matcher.set_hard_threshold({min_lambda, max_lambda});
    matcher.set_candidates(std::move(candidate_subsets),
                           std::vector<double>(statistics.begin(), statistics.end()));
    matcher.set_interrupts_checker([]() { Rcpp::checkUserInterrupt(); });

    auto counts = matrix_to_counts(*case_counts);
    case_counts.reset();
    auto result = matcher.match(counts, min_controls, min_call_rate);

    List ret;
    NumericVector lambda(result.lambdas.begin(), result.lambdas.end());
    NumericVector pvals(result.pvals.begin(), result.pvals.end());
    NumericVector stats(result.statistics.begin(), result.statistics.end());
    IntegerVector names(result.lambda_i.begin(), result.lambda_i.end());
    IntegerVector pvals_num(result.pvals_num.begin(), result.pvals_num.end());
    IntegerVector optimal_clusters(result.optimal_prefix.begin(), result.optimal_prefix.end());
    NumericVector optimal_lambda = {result.optimal_lambda};

    lambda.attr("names") = names;
    stats.attr("names") = names;
    pvals_num.attr("names") = names;
    ret["lambda"] = lambda;
    ret["optimal_lambda"] = optimal_lambda;
    ret["statistics"] = stats;
    ret["controls"] = optimal_clusters + 1;
    ret["pvals"] = pvals;
    ret["snps"] = pvals_num;
    return ret;
}
