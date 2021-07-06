#include "include/subsample.h"

#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include "include/third-party/cxxpool.h"

namespace mvn {

subsample::subsample(std::shared_ptr<const mvn::Matrix> X, const Clustering& clst, const mvn::Vector& mean,
                     const mvn::Matrix& cov)
        : test{std::make_shared<mvn_test>(mvn_test(X, clst, cov, mean))},
          clst(clst),
          wheel(std::random_device()()) {
}

void subsample::run(size_t iterations, size_t restarts, double t_start, double c, size_t pool_size, size_t start, size_t size_ub,
                    size_t step) {
    using std::vector;
    using namespace std::chrono_literals;

    best.clear();
    best_scores.clear();

    cxxpool::thread_pool pool(pool_size);
    size_t curr_size = start;

    while (curr_size <= size_ub) {
        Rcpp::checkUserInterrupt();
        std::vector<std::future<std::shared_ptr<mvn_test>>> thread_solutions;
        for (size_t t = 0; t < restarts; t++) {
            thread_solutions.push_back(pool.push([test = std::as_const(test), iterations, t_start, curr_size, c,
                                                         seed = wheel()]() -> std::shared_ptr<mvn_test> {
                double t = t_start;
                std::mt19937 mersenne_wheel(seed);
                std::shared_ptr<mvn_test> local_test = test->clone();
                while (local_test->subsample_size() < curr_size) {
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

        for (size_t i = 0; i < thread_solutions.size(); i++) {
            auto& future = thread_solutions[i];
            std::shared_ptr<mvn_test> thread = future.get();
            if (i == 0) {
                test = thread;
            } else if (*thread < *test) {
                test = thread;
            }
        }

        best.push_back(test->current_subset());
        best_scores.push_back(test->get_normality_statistic());

        curr_size += step;

        std::ofstream fout("sets/" + std::to_string(test->subsample_size()) + ".txt");
        for (int gr: test->current_subset()) {
            for (int el: clst.elements(gr)) {
                fout << el << std::endl;
            }
        }

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

}