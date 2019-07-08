#include "include/utils.h"
#include <Rcpp.h>
#include <algorithm>

using std::vector;
using std::tuple;

const double chi2_bdry = 10.82757;

Clustering::Iterator::Iterator(const Clustering& object) :obj(object) {}

Control Clustering::Iterator::next() {
    int curr_size = obj.cluster_sizes[obj.permutation[cluster]];
    Control ret;
    ret.number = obj.clusters[obj.permutation[cluster]][shift];
    ++shift;
    if (shift == curr_size) {
        ret.last = true;
        cluster++;
        shift = 0;
    }
    return ret;
}


Clustering::Clustering(const vector<double>& residuals, const vector<int>& clustering) {
    n = residuals.size();
    if (clustering.empty()) {
        throw std::invalid_argument("Clustering must not be empty");
    }

    int n_clsuters = *std::max_element(clustering.begin(), clustering.end()) + 1;
    cluster_sizes.resize(n_clsuters);
    res.resize(n_clsuters);
    clusters.resize(n_clsuters);
    for (int i = 0; i < clustering.size(); i++) {
        int cl = clustering[i];
        ++cluster_sizes[cl];
        res[cl] += residuals[i];
        clusters[cl].push_back(i);
    }
    for (int i = 0; i < n_clsuters; i++) {
        res[i] /= cluster_sizes[i];
    }

    permutation.resize(n_clsuters);
    std::iota(permutation.begin(), permutation.end(), 0);
    std::sort(permutation.begin(), permutation.end(), [this](int pos_a, int pos_b) {
        return res[pos_a] < res[pos_b];
    });
}

Clustering::Iterator Clustering::iterator() const {
    return Clustering::Iterator(*this);
}

size_t Clustering::size() const {
    return n;
}

std::vector<bool> check_user_counts(std::vector<std::vector<int>>& case_counts) {
    std::vector<bool> mask;
    for (auto& case_count : case_counts) {
        int homref = case_count[0];
        int het = case_count[1];
        int hom = case_count[2];
        mask.push_back(check_counts(homref, het, hom));
    }
    return mask;
}

double chi2_aux(double obs, double exp) {
    double chi = std::abs(obs - exp) - 0.5;
    return (chi * chi) / exp;
}

double chi2(int k, int n, std::function<double(double)>& qchisq) {
    double a = n <= 10 ? 3.0 / 8.0 : 0.5;
    return qchisq((k + 1 - a) / (n + 1 - 2 * a));
}

matching_results::matching_results(int opt_prefix, double opt_lmd, std::vector<double>&& p_values,
        std::vector<double>&& lmbds, std::vector<int>&& lmbd_i, std::vector<int>&& pvals_number)
        :optimal_prefix(opt_prefix), optimal_lambda(opt_lmd), pvals(std::move(p_values)),
        lambdas(std::move(lmbds)), lambda_i(std::move(lmbd_i)), pvals_num(std::move(pvals_number)) {}

double lambda_1000(double lambda_observed, long cases, long controls) {
    double coefficient = 1.0 / cases + 1.0 / controls;
    return 1.0 + ((lambda_observed - 1) * coefficient) / (2 * 0.001);
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

matching_results select_controls_impl(vector<vector<int>>& gmatrix, const Clustering& clustering,
                                                  vector<vector<int>>& case_counts,
                                                  std::function<double(double)>& qchisq,
                                                  double min_lambda, double lb_lambda,
                                                  double max_lambda, double ub_lambda,
                                                  size_t min_controls) {
    min_controls = std::max(min_controls, 1UL);
    std::vector<bool> snp_mask = check_user_counts(case_counts);

    auto it = clustering.iterator();
    unsigned long n = clustering.size();
    unsigned long m = case_counts.size();
    std::vector<std::vector<int>> counts(m, vector<int>(3));
    std::vector<lm> lms;
    unsigned rank = 0;
    for (size_t i = 0; i < m; i++) {
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

    for (size_t k = 0; k < n; k++) {
        auto instance = it.next();
        int i = instance.number;
        if (i % 100 == 0) {
            Rcpp::checkUserInterrupt();
        }

        ++rank;
        std::vector<double> pvals;
        for (size_t j = 0; j < m; j++) {
            if (!snp_mask[j]) {
                continue;
            }
            int cur = gmatrix[i][j];
            std::vector<int>& cts = counts[j];
            if (cur != -1) {
                ++cts[cur];
                lms[j].set(cur, cur, 0, cts[cur]);
            }
            if (k >= min_controls - 1 && instance.last) {
                if (!check_counts(cts[0], cts[1], cts[2])) {
                    continue;
                }
                lms[j].solve();
                pvals.push_back(pval_t(lms[j].compute_t(rank - 2), rank - 2));
            }
        }
        unsigned long n_pvals = pvals.size();
        if (k >= min_controls - 1 && pvals.size() > 1 && instance.last) {
            lm pvals_lm(n_pvals);
            std::sort(pvals.begin(), pvals.end());
            for (size_t j = 0; j < n_pvals; j++) {
                pvals_lm.set(j, chi2(j, n_pvals, qchisq), qchisq(pvals[j]), 1, false);
            }
            pvals_lm.solve();

            double cur_lambda = pvals_lm.get_lambda();
            long cases = std::accumulate(case_counts[0].begin(), case_counts[0].end(), 0);
            long controls = std::accumulate(counts[0].begin(), counts[0].end(), 0);
            cur_lambda = lambda_1000(cur_lambda, cases, controls);

            lambdas.push_back(cur_lambda);
            lambda_i.push_back(k + 1);
            pvals_num.push_back(pvals.size());
            double lambda_dist = std::max(lambda - ub_lambda, lb_lambda - lambda);
            double cur_lambda_dist = std::max(cur_lambda - ub_lambda, lb_lambda - cur_lambda);
            if (cur_lambda < max_lambda && cur_lambda > min_lambda) {
                if ((cur_lambda < ub_lambda && cur_lambda > lb_lambda) || cur_lambda_dist < lambda_dist) {
                    optimal_prefix = k + 1;
                    lambda = cur_lambda;
                    optimal_pvals = pvals;
                }
            }
        }
    }

    return matching_results(optimal_prefix, lambda, std::move(optimal_pvals), std::move(lambdas), std::move(lambda_i),
                            std::move(pvals_num));
}

