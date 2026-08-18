#ifndef SRC_MATCHING_H
#define SRC_MATCHING_H

#include <cstdint>
#include <vector>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <functional>
#include <stdexcept>

#include "mvn_test.h"
#include "lm.h"


namespace matching {

class Counts {
    int counts[3];
public:
    Counts() :counts{0, 0, 0} {}
    Counts(int hom, int het, int alt) :counts{hom, het, alt} {}

    const int& operator[] (size_t i) const {
        check_bounds(i);
        return counts[i];
    }
    int& operator[] (size_t i) {
        check_bounds(i);
        return counts[i];
    }

    int sum() const {
        return counts[0] + counts[1] + counts[2];
    }

private:
    void check_bounds(size_t i) const {
        if (i > 2) {
            throw std::domain_error("Out of bound");
        }
    }
};

class ClusterCounts {
    std::uint8_t counts[3];
public:
    ClusterCounts() :counts{0, 0, 0} {}

    const std::uint8_t& operator[] (size_t i) const {
        check_bounds(i);
        return counts[i];
    }
    std::uint8_t& operator[] (size_t i) {
        check_bounds(i);
        return counts[i];
    }

private:
    void check_bounds(size_t i) const {
        if (i > 2) {
            throw std::domain_error("Out of bound");
        }
    }
};

class lambda_range {
    double lb;
    double ub;
public:
    lambda_range(double lb, double ub);
    bool in(double lambda);
    lambda_range();
    double distance(double lambda);
};

struct matching_results {
    std::vector<int> optimal_prefix;
    std::vector<double> pvals;
    std::vector<double> lambdas;
    std::vector<double> statistics;
    std::vector<int> lambda_i;
    std::vector<int> pvals_num;
    double optimal_lambda;

    matching_results(std::vector<int>&& optimal_prefix, std::vector<double>&& p_values, std::vector<double>&& lmbds,
                     std::vector<double>&& statistics, std::vector<int>&& lmbd_i, std::vector<int>&& pvals_number,
                     double oprimal_lambda);
};

// Scores candidate control subsets genetically: per-variant counts, call rate,
// association p-values and the resulting lambda_GC. The candidates themselves come
// from the multivariate normality search, which knows nothing about genetics and
// runs separately (see mvn::subsample).
class matching {
    static constexpr double EPS = 1e-6;

    std::vector<std::vector<ClusterCounts>> cluster_counts;

    mvn::Clustering clustering;

    // One candidate subset of clusters per target size, each with the normality
    // statistic the search reached for it.
    std::vector<std::vector<size_t>> candidates;
    std::vector<double> candidate_statistics;

    std::function<void()> interrupts_checker;
    std::function<double(double)> qchisq;

    lambda_range hard_threshold;
    lambda_range soft_threshold;
public:
    matching(std::vector<std::vector<ClusterCounts>>&& cluster_counts, mvn::Clustering clustering);
    void set_candidates(std::vector<std::vector<size_t>>&& candidates, std::vector<double>&& statistics);
    void set_qchi_sq_function(const std::function<double(double)>& f);
    matching_results match(const std::vector<Counts>& case_counts, unsigned min_controls = 1, double min_call_rate = 0.95);

    void set_interrupts_checker(std::function<void()> checker) {
        interrupts_checker = checker;
    }

    void set_soft_threshold(lambda_range range);
    void set_hard_threshold(lambda_range range);
private:
    double get_lambda(std::vector<double>& pvals);

    Counts count_controls(const std::vector<size_t>& groups, size_t variant);
};

}

#endif //SRC_MATCHING_H
