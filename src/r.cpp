#include <Rcpp.h>
#include <vector>

#include "include/qchisq.h"
#include "include/utils.h"

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
        ret[i] = check_counts(homref, het, hom, af, ac, boundary);
    }
    return ret;
}

// [[Rcpp::export]]
List select_controls_cpp(IntegerMatrix& gmatrix, NumericVector& residuals,
                     IntegerMatrix& cc, IntegerVector& clustering,
                     NumericVector& chi2fn,
                     NumericVector min_lambda, NumericVector lb_lambda,
                     NumericVector max_lambda, NumericVector ub_lambda, IntegerVector min) {
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

    vector<int> clust_vec(clustering.begin(), clustering.end());

    for (int i = 0; i < n_c; i++) {
        resid[i] = residuals[i];
        for (int j = 0; j < gmatrix.nrow(); j++) {
            int value = gmatrix(j, i);
            gmatrix_c[i][j] = value == NA_INTEGER ? -1 : value;
        }
    }
    for (int i = 0; i < n_snps; i++) {
        for (int j = 0; j < 3; j++) {
            case_counts[i][j] = cc(i, j);
        }
    }

    int min_controls = min[0];
    std::function<double(double)> qchi = q.function();
    Clustering cl(resid, clust_vec);
    auto result = select_controls_impl(gmatrix_c, cl, case_counts, qchi, min_lambda[0], lb_lambda[0],
            max_lambda[0], ub_lambda[0], min_controls);

    std::vector<int> permutation;
    auto it = cl.iterator();
    for (int i = 0; i < cl.size(); i++) {
        permutation.push_back(it.next().number);
    }

    List ret;
    NumericVector lambda(result.lambdas.begin(), result.lambdas.end());
    IntegerVector n_controls(1, result.optimal_prefix);
    NumericVector pvals(result.pvals.begin(), result.pvals.end());
    NumericVector optimal_lambda(1, result.optimal_lambda);
    IntegerVector names(result.lambda_i.begin(), result.lambda_i.end());
    IntegerVector pvals_num(result.pvals_num.begin(), result.pvals_num.end());
    IntegerVector permutationR(permutation.begin(), permutation.end());

    lambda.attr("names") = names;
    pvals_num.attr("names") = names;
    ret["lambda"] = lambda;
    ret["controls"] = n_controls;
    ret["pvals"] = pvals;
    ret["optimal_lambda"] = optimal_lambda;
    ret["snps"] = pvals_num;
    ret["permutation"] = permutationR;
    return ret;
}
