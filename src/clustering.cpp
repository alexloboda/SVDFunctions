#include "include/sskm.h"
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::IntegerVector sskm_cpp(const Rcpp::NumericMatrix& X, const Rcpp::IntegerVector k, const Rcpp::IntegerVector max_iter, 
                             const Rcpp::NumericVector tol) {
    clustering::SameSizeKMeans sskm(k[0], max_iter[0], tol[0]);
    std::vector<std::vector<double>> X_vec(X.nrow());
    for (int i = 0; i < X.nrow(); i++) {
        X_vec[i].resize(X.ncol());
        for (int j = 0; j < X.ncol(); j++) {
            X_vec[i][j] = X(i, j);
        }
    }
    sskm.fit(X_vec);
    return Rcpp::wrap(sskm.get_labels());
}