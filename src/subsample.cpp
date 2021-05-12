#include "include/subsample.h"

#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include "include/third-party/cxxpool.h"

namespace mvn {

subsample::subsample(std::shared_ptr<const mvn::Matrix> X, const Clustering& clst, const mvn::Vector& mean, const mvn::Matrix& cov)
        :test{std::make_shared<mvn_test_fixed>(mvn_test_fixed(X, clst, cov, mean))},
         clst(clst),
         wheel(std::random_device()()) {
    Rcpp::Rcerr << "Test is prepared, mean vector is " << mean << std::endl <<
                "Covariance matrix is :" << cov << std::endl;
    Rcpp::Rcerr << std::flush;
}



void subsample::run(size_t iterations, size_t restarts, double t_start, double c, int pool_size, int start, int size_ub, int step) {
    using std::vector;
    using namespace std::chrono_literals;

    best.clear();
    best_scores.clear();

    cxxpool::thread_pool pool(pool_size);
    int curr_size = start;

    while (curr_size <= size_ub) {
        Rcpp::checkUserInterrupt();
        std::vector<std::future<std::shared_ptr<mvn_test>>> thread_solutions;
        for (int t = 0; t < restarts; t++) {
            thread_solutions.push_back(pool.push([test = std::as_const(test), iterations, t_start, curr_size, c,
                                                         seed = wheel()]() -> std::shared_ptr<mvn_test> {
                double t = t_start;
                std::mt19937 mersenne_wheel(seed);
                std::shared_ptr<mvn_test> local_test = test->clone();
                while (test->subsample_size() < curr_size) {
                    local_test->add_one();
                }
                std::uniform_real_distribution<double> random_unif(0.0, 1.0);
                double score = local_test->get_normality_statistic();
                for (size_t k = 0; k < iterations; k++) {
                    t = c * t;
                    local_test->swap_once();
                    double new_score = local_test->get_normality_statistic();
                    if (new_score < score) {
                        score = new_score;
                    } else {
                        double p = std::exp((score - new_score) / t);
                        if (random_unif(mersenne_wheel) >= p) {
                            local_test->swap_once(true);
                        } else {
                            score = new_score;
                        }
                    }
                }
                return local_test;
            }));
        }

        for (auto& future: thread_solutions) {
            std::shared_ptr<mvn_test> thread = future.get();
            if (*thread < *test) {
                test = thread;
            }
        }

        best.push_back(test->current_subset());
        best_scores.push_back(test->get_normality_statistic());

        curr_size += step;

        if (curr_size > test->sample_size()) {
            return;
        }
    }
}

size_t subsample::solutions() const {
    return best.size();
}

std::vector<size_t> subsample::get_solution(size_t k) const {
    return best.at(k);
}

subsample::subsample() {}

subsample::subsample(std::shared_ptr<const Matrix> X, const Clustering& clst)
    :test{std::make_shared<mvn_test_gen>(mvn_test_gen(X, clst))},
    clst(clst),
    wheel(std::random_device()()) {}

}

