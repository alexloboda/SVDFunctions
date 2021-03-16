#ifndef SRC_MATCHING_H
#define SRC_MATCHING_H

#include <vector>
#include <eigen3/Eigen/Dense>
#include <functional>
#include <stdexcept>

#include "subsample.h"
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
    std::vector<int> lambda_i;
    std::vector<int> pvals_num;

    matching_results(std::vector<int>&& optimal_prefix, std::vector<double>&& p_values, std::vector<double>&& lmbds,
                     std::vector<int>&& lmbd_i, std::vector<int>&& pvals_number);
};

class matching {
    static constexpr double EPS = 1e-6;

    Eigen::MatrixXi controls_gmatrix;
    Eigen::MatrixXd controls_space;

    mvn::Clustering clustering;
    mvn::subsample subsampling;

    std::function<void()> interrupts_checker;
    std::function<double(double)> qchisq;

    lambda_range hard_threshold;
    lambda_range soft_threshold;
public:
    matching(Eigen::MatrixXi controls, Eigen::MatrixXd space, mvn::Clustering clustering);
    void process_mvn(const Eigen::MatrixXd& directions, Eigen::VectorXd mean, int threads);
    void set_qchi_sq_function(const std::function<double(double)>& f);
    matching_results match(const std::vector<Counts>& case_counts, unsigned min_controls = 1);

    void set_interrupts_checker(std::function<void()> checker) {
        interrupts_checker = checker;
    }

    void set_soft_threshold(lambda_range range);
    void set_hard_threshold(lambda_range range);
private:
    Eigen::VectorXd cluster_mean(int cluster) const;
    double get_lambda(std::vector<double>& pvals);

    Counts count_controls(const std::vector<int>& vector, size_t j);
};

}

#endif //SRC_MATCHING_H
