#include <Rcpp.h>
#include <vector>
#include <eigen3/Eigen/Dense>
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

// [[Rcpp::export]]
List subsample_mvn(NumericMatrix& matrix, IntegerVector size) {
    //std::vector<int> clusters(matrix.ncol());
    //std::iota(clusters.begin(), clusters.end(), 0);
    //mvn::Clustering clustering(clusters);
    //mvn::subsample annealing(r_to_cpp(matrix), clustering);
    //annealing.run(10'000, 4, 1.0, 0.995, std::thread::hardware_concurrency(), size[0], size[0], 1);

    //List ret;
    //auto points = annealing.get_solution(0);
    //IntegerVector subset(points.begin(), points.end());
    //ret["points"] = subset + 1;

    //return ret;
}

// [[Rcpp::export]]
List select_controls_cpp(IntegerMatrix& gmatrix,
                     NumericMatrix& gmatrix_rs,
                     NumericVector& mean, NumericMatrix& directions,
                     IntegerMatrix& cc, IntegerVector& clustering,
                     NumericVector& chi2fn,
                     NumericVector min_lambda, NumericVector lb_lambda,
                     NumericVector max_lambda, NumericVector ub_lambda,
                     IntegerVector min, IntegerVector max, IntegerVector step) {
    vector<double> precomputed_chi(chi2fn.begin(), chi2fn.end());
    qchi2 q(precomputed_chi);

    auto gmatrix_counts = r_to_cpp(gmatrix);
    auto case_counts = r_to_cpp(cc);
    auto principal_directions = r_to_cpp(directions);
    auto gm_rs = r_to_cpp(gmatrix_rs);

    vector<int> clust_vec(clustering.begin(), clustering.end());

    int min_controls = min[0];
    int max_controls = max[0];
    int step_clusters = step[0];
    mvn::Clustering cl(clust_vec);

    matching::matching matcher(gmatrix_counts, gm_rs, cl);
    matcher.set_qchi_sq_function(q.function());
    matcher.set_soft_threshold({lb_lambda[0], ub_lambda[0]});
    matcher.set_hard_threshold({min_lambda[0], max_lambda[0]});
    matcher.process_mvn(*principal_directions, r_to_cpp(mean), std::thread::hardware_concurrency(),
                        min_controls, max_controls, step_clusters);
    matcher.set_interrupts_checker([]() { Rcpp::checkUserInterrupt(); });

    auto result = matcher.match(matrix_to_counts(*case_counts), min_controls);

    List ret;
    NumericVector lambda(result.lambdas.begin(), result.lambdas.end());
    NumericVector pvals(result.pvals.begin(), result.pvals.end());
    IntegerVector names(result.lambda_i.begin(), result.lambda_i.end());
    IntegerVector pvals_num(result.pvals_num.begin(), result.pvals_num.end());
    IntegerVector optimal_controls(result.optimal_prefix.begin(), result.optimal_prefix.end());

    lambda.attr("names") = names;
    pvals_num.attr("names") = names;
    ret["lambda"] = lambda;
    ret["controls"] = optimal_controls + 1;
    ret["pvals"] = pvals;
    ret["snps"] = pvals_num;
    return ret;
}
