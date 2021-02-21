#include "include/subsample.h"

#include <Rcpp.h>
#include <numeric>
#include <iostream>
#include <random>

namespace mvn {

subsample::subsample(const mvn::Matrix& X, const mvn::Vector& mean, const mvn::Matrix& cov)
        :test(X, cov, mean),
         wheel(std::random_device()()),
         random_unif(0.0, 1.0) {}

std::vector<size_t> subsample::run(size_t iterations, size_t restarts, double t0) {
    using std::vector;

    best.clear();
    best_scores.clear();

    while (test.subsample_size() < test.sample_size()) {
        double score = test.get_normality_statistic();

        best.push_back(test.current_subset());
        best_scores.push_back(score);

        for (size_t i = 0; i < restarts; i++) {
            Rcpp::checkUserInterrupt();
            for (size_t k = 1; k <= iterations; k++) {
                double t = t0 / (1.0 + std::log((double) k));

                test.swap_once();
                double new_score = test.get_normality_statistic();
                if (new_score < score) {
                    score = new_score;
                    if (score < best_scores.back()) {
                        best_scores.back() = score;
                        best.back() = test.current_subset();
                    }
                } else {
                    double p = std::exp((score - new_score) / t);
                    if (random_unif(wheel) >= p) {
                        test.swap_once(true);
                    } else {
                        score = new_score;
                    }
                }
            }
            t0 /= 2.0;
        }
        test.add_one();
    }
}

size_t subsample::min_size() const {
    return best.front().size();
}

size_t subsample::max_size() const {
    return best.back().size();
}

double subsample::get_statistic(size_t size) const {
    return best_scores.at(size - min_size());
}

std::vector<size_t> subsample::get_solution(size_t size) const {
    return best.at(size - min_size());
}

subsample::subsample() {}

}

