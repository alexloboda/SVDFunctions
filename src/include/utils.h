#include <utility>

#ifndef SRC_UTILS_H
#define SRC_UTILS_H

#include <vector>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <functional>

#include "lm.h"

using std::vector;
using std::tuple;

const double chi2_bdry = 10.82757;

namespace {

void sort_controls(vector<vector<int>>& gmatrix, vector<double>& residuals) {
    vector<int> perm(residuals.size());
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(), [&residuals](int pos_a, int pos_b){
        return residuals[pos_a] < residuals[pos_b];
    });
    vector<vector<int>> gm;
    vector<double> sorted_residuals;
    for (int i = 0; i < residuals.size(); i++) {
        gm.push_back(std::move(gmatrix[perm[i]]));
        sorted_residuals.push_back(residuals[perm[i]]);
    }
    gmatrix = std::move(gm);
    residuals = std::move(sorted_residuals);
}

double chi2_aux(double obs, double exp) {
    double chi = std::abs(obs - exp) - 0.5;
    return (chi * chi) / exp;
}

bool check_counts(unsigned hom_ref, unsigned het, unsigned hom) {
    if (hom < hom_ref) {
        std::swap(hom, hom_ref);
    }
    double AC = 2 * hom_ref + het;
    unsigned n = hom_ref + het + hom;
    double p = AC / (2 * n);
    double q = 1 - p;
    double chi2 = chi2_aux(hom_ref, n * p * p) + chi2_aux(het, 2 * n * p * q) + chi2_aux(hom, n * q * q);
    return p > 0.05 &&  AC > 10 && chi2 < chi2_bdry;
}

std::vector<bool> check_user_counts(vector<vector<int>>& case_counts) {
    std::vector<bool> mask;
    for (int i = 0; i < case_counts.size(); i++) {
        int homref = case_counts[i][0];
        int het = case_counts[i][1];
        int hom = case_counts[i][2];
        mask.push_back(check_counts(homref, het, hom));
    }
    return mask;
}

double chi2(int k, int n, std::function<double(double)>& qchisq) {
    double a = n <= 10 ? 3.0 / 8.0 : 0.5;
    return qchisq((k + 1 - a) / (n + 1 - 2 * a));
}

}

struct matching_results {
    int optimal_prefix;
    double optimal_lambda;
    vector<double> pvals;
    vector<double> lambdas;
    vector<int> lambda_i;
    vector<int> pvals_num;

    matching_results(int opt_prefix, double opt_lmd, vector<double>&& p_values, vector<double>&& lmbds,
                     vector<int>&& lmbd_i, vector<int>&& pvals_number)
            :optimal_prefix(opt_prefix), optimal_lambda(opt_lmd), pvals(std::move(p_values)), 
	    lambdas(std::move(lmbds)), lambda_i(std::move(lmbd_i)), pvals_num(std::move(pvals_number)) {}
};

matching_results select_controls_impl(vector<vector<int>>& gmatrix, vector<double>& residuals,
                                                  vector<vector<int>>& case_counts, std::function<double(double)>& qchisq,
                                                  double min_lambda, double lb_lambda,
                                                  double max_lambda, double ub_lambda,
                                                  int min_controls = 500, int bin = 1) {
    bin = std::max(bin, 1);
    std::vector<bool> snp_mask = check_user_counts(case_counts);

    sort_controls(gmatrix, residuals);
    unsigned long n = residuals.size();
    unsigned long m = case_counts.size();
    std::vector<std::vector<int>> counts(m, vector<int>(3));
    std::vector<lm> lms;
    unsigned rank = 0;
    for (int i = 0; i < m; i++) {
        lm model(6);
        for (int j = 0; j < 3; j++) {
            model.set(3 + j, j, 1, case_counts[i][j]);
            if (i == 0) {
                rank += case_counts[i][j];
            }
        }
        lms.push_back(model);
    }

    double lambda = std::numeric_limits<double>::infinity();
    std::vector<double> optimal_pvals;
    std::vector<double> lambdas;
    std::vector<int> lambda_i;
    std::vector<int> pvals_num;
    int optimal_prefix = -1;

    for (int i = 0; i < n; i++) {
        if (i % 100 == 0) {
            Rcpp::checkUserInterrupt();
        }

        ++rank;
        std::vector<double> pvals;
        for (int j = 0; j < m; j++) {
            if (!snp_mask[j]) {
                continue;
            }
            int cur = gmatrix[i][j];
            std::vector<int>& cts = counts[j];
            ++cts[cur];
            lms[j].set(cur, cur, 0, cts[cur]);
            if (i >= min_controls - 1 && (i + 1) % bin == 0) {
                if (!check_counts(cts[0], cts[1], cts[2])) {
                    continue;
                }
                lms[j].solve();
                pvals.push_back(pval_t(lms[j].compute_t(rank - 2), rank - 2));
            }
        }
        unsigned long n_pvals = pvals.size();
        if (i >= min_controls - 1 && !pvals.empty() && (i + 1) % bin == 0) {
            lm pvals_lm(n_pvals);
            std::sort(pvals.begin(), pvals.end());
            for (int j = 0; j < n_pvals; j++) {
                pvals_lm.set(j, chi2(j, n_pvals, qchisq), qchisq(pvals[j]), 1, false);
            }
            pvals_lm.solve();
            double cur_lambda = pvals_lm.get_lambda();
            lambdas.push_back(cur_lambda);
            lambda_i.push_back(i + 1);
            pvals_num.push_back(pvals.size());
            double lambda_dist = std::max(lambda - ub_lambda, lb_lambda - lambda);
            double cur_lambda_dist = std::max(cur_lambda - ub_lambda, lb_lambda - cur_lambda);
            if (cur_lambda < max_lambda && cur_lambda > min_lambda) {
                if ((cur_lambda < ub_lambda && cur_lambda > lb_lambda) || cur_lambda_dist < lambda_dist) {
                    optimal_prefix = i + 1;
                    lambda = cur_lambda;
                    optimal_pvals = pvals;
                }
            }
        }
    }

    return matching_results(optimal_prefix, lambda, std::move(optimal_pvals), std::move(lambdas), std::move(lambda_i),
                            std::move(pvals_num));
}

#endif
