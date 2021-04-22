#include "include/subsample.h"

#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include "include/third-party/cxxpool.h"

namespace mvn {

subsample::subsample(const mvn::Matrix& X, const Clustering& clst, const mvn::Vector& mean, const mvn::Matrix& cov)
        :test(X, clst, cov, mean),
         clst(clst),
         wheel(std::random_device()()) {
    Rcpp::Rcerr << "Test is prepared, mean vector is " << mean << std::endl <<
                "Covariance matrix is :" << cov << std::endl;
    Rcpp::Rcerr << std::flush;
}

void subsample::run(size_t iterations, size_t restarts, double t_start, int pool_size, int start, int size_ub, int step) {
    using std::vector;
    using namespace std::chrono_literals;

    best.clear() ;
    best_scores.clear();
    pool_size = 2;
    cxxpool::thread_pool pool(pool_size);

    while (test.subsample_size() < start) {
        test.add_one();
    }

    while (test.subsample_size() <= size_ub) {
        double t0 = t_start;
        for (size_t i = 0; i < restarts; i++) {
            Rcpp::checkUserInterrupt();
            std::vector<std::future<mvn_test>> thread_solutions;
            for (int t = 0; t < pool_size; t++) {
                thread_solutions.push_back(pool.push([&test = std::as_const(test), iterations, t0,
                                                             seed = wheel()]() -> mvn_test {
                    std::mt19937 mersenne_wheel(seed);
                    std::uniform_real_distribution<double> random_unif(0.0, 1.0);
                    mvn_test local_test = test;
                    double score = local_test.get_normality_statistic();
                    for (size_t k = 1; k <= iterations; k++) {
                        double t = t0 / (1.0 + std::log((double) k));
                        local_test.swap_once();
                        double new_score = local_test.get_normality_statistic();
                        if (new_score < score) {
                            score = new_score;
                        } else {
                            double p = std::exp((score - new_score) / t);
                            if (random_unif(mersenne_wheel) >= p) {
                                local_test.swap_once(true);
                            } else {
                                score = new_score;
                            }
                        }
                    }
                    return local_test;
                }));
            }

            for (auto& future: thread_solutions) {
                test = std::min(test, future.get());
            }

            t0 /= 2.0;
        }

        best.push_back(test.current_subset());
        best_scores.push_back(test.get_normality_statistic());

        Rcpp::Rcerr << "Subsampling set of size " << test.subsample_size();
        Rcpp::Rcerr << ". Current score is " << test.get_normality_statistic() << std::endl;

        std::ofstream fout("sets/" + std::to_string(test.subsample_size()) + ".txt");
        for (int gr: test.current_subset()) {
            for (int el: clst.elements(gr)) {
                fout << el << std::endl;
            }
        }

        if (test.subsample_size() == test.sample_size()) {
            return;
        }
        for (int i = 0; i < step; i++) {
            test.add_one();
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

