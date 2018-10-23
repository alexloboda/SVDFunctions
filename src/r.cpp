#include <Rcpp.h>
#include <vector>

#include "qchisq.h"
#include "utils.h"

using namespace Rcpp;
using std::vector;

// [[Rcpp::export]]
List select_controls_cpp(NumericMatrix& gmatrix, NumericVector& residuals,
                     NumericMatrix& cc, NumericVector& chi2fn, NumericVector min,
                     IntegerVector bin_size) {
    int n_c = gmatrix.ncol();
    int n_snps = gmatrix.nrow();
    vector<double> precomputed_chi(chi2fn.length());
    for (int i = 0; i < precomputed_chi.size(); i++) {
        precomputed_chi[i] = chi2fn[i];
    }
    qchi2 q(precomputed_chi);

    vector<vector<int>> gmatrix_c(n_c, vector<int>(n_snps));
    vector<vector<int>> case_counts(n_snps, vector<int>(3));
    vector<double> resid(n_c);

    for (int i = 0; i < n_c; i++) {
        resid[i] = residuals[i];
        for (int j = 0; j < gmatrix.nrow(); j++) {
            gmatrix_c[i][j] = (int)std::lround(gmatrix(j, i));
        }
    }
    for (int i = 0; i < n_snps; i++) {
        for (int j = 0; j < 3; j++) {
            case_counts[i][j] = (int)std::lround(cc(i, j));
        }
    }
    int min_controls = (int)std::lround(min[0]);
    int bin = (int)std::lround(bin_size[0]);
    std::function<double(double)> qchi = q.function();
    auto result = select_controls_impl(gmatrix_c, resid, case_counts, qchi, min_controls, bin);
    List ret;
    NumericVector lambda(result.lambdas.begin(), result.lambdas.end());
    IntegerVector n_controls(1, result.optimal_prefix);
    NumericVector pvals(result.pvals.begin(), result.pvals.end());
    NumericVector optimal_lambda(1, result.optimal_lambda);
    IntegerVector names(result.lambda_i.begin(), result.lambda_i.end());
    lambda.attr("names") = names;
    ret["lambda"] = lambda;
    ret["controls"] = n_controls;
    ret["pvals"] = pvals;
    ret["optimal_lambda"] = optimal_lambda;
    return ret;
}
