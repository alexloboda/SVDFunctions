#include <utility>

#ifndef SRC_UTILS_H
#define SRC_UTILS_H

#include <vector>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <functional>

#include "lm.h"

bool check_counts(unsigned hom_ref, unsigned het, unsigned hom);

struct matching_results {
    int optimal_prefix;
    double optimal_lambda;
    std::vector<double> pvals;
    std::vector<double> lambdas;
    std::vector<int> lambda_i;
    std::vector<int> pvals_num;

    matching_results(int opt_prefix, double opt_lmd, std::vector<double>&& p_values, std::vector<double>&& lmbds,
                     std::vector<int>&& lmbd_i, std::vector<int>&& pvals_number);
};

struct Control {
    int number;
    bool last = false;
};

class Clustering {
    std::vector<int> cluster_sizes;
    std::vector<double> res;
    std::vector<std::vector<int>> clusters;
    std::vector<int> permutation;
    size_t n;
public:
    class Iterator {
        const Clustering &obj;
        int cluster = 0;
        int shift = 0;
    public:
        explicit Iterator(const Clustering& object);
        Control next();
    };
    Clustering(const std::vector<double>& residuals, const std::vector<int>& clustering);
    Iterator iterator() const;
    size_t size() const;
};

matching_results select_controls_impl(std::vector<std::vector<int>>& gmatrix, const Clustering& clustering,
                                                  std::vector<std::vector<int>>& case_counts,
                                                  std::function<double(double)>& qchisq, double min_lambda,
                                                  double lb_lambda, double max_lambda, double ub_lambda,
                                                  std::size_t min_controls = 500);

#endif
