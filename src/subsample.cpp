#include "include/subsample.h"

#include <Rcpp.h>
#include <iostream>
#include <random>
#include "include/third-party/cxxpool.h"

namespace {
    class solution {
        std::vector<size_t> elements;
        double score;
    public:
        solution(std::vector<size_t> elements, double score) :elements(elements), score(score) {}
        solution& operator=(const solution&) = default;
        solution() :score(std::numeric_limits<double>::infinity()){}
        bool operator<(const solution& other) const {
            return score < other.score;
        }
        std::vector<size_t> elems() const {
            return elements;
        }
        double fitness() const {
            return score;
        }
    };

}

namespace mvn {

subsample::subsample(const mvn::Matrix& X, const Clustering& clst, const mvn::Vector& mean, const mvn::Matrix& cov)
        :test(X, clst, cov, mean),
         wheel(std::random_device()()) {}

void subsample::run(size_t iterations, size_t restarts, double t0, int pool_size) {
    using std::vector;

    best.clear() ;
    best_scores.clear();
    cxxpool::thread_pool pool(pool_size);

    while (test.subsample_size() < test.sample_size()) {
        Rcpp::Rcerr << "Subsampling set of size " << test.subsample_size();
        Rcpp::Rcerr << ". Current score is " << test.get_normality_statistic() << std::endl;

        for (size_t i = 0; i < restarts; i++) {
            Rcpp::checkUserInterrupt();
            std::vector<std::future<solution>> thread_solutions;
            for (size_t t = 0; t < pool_size; t++) {
                thread_solutions.push_back(pool.push([&test = std::as_const(test), iterations, t0,
                                                             seed = wheel()]() -> solution {
                    std::mt19937 mersenne_wheel(seed);
                    std::uniform_real_distribution<double> random_unif(0.0, 1.0);
                    solution best_thread;
                    mvn_test local_test = test;
                    double score = local_test.get_normality_statistic();
                    for (size_t k = 1; k <= iterations; k++) {
                        double t = t0 / (1.0 + std::log((double) k));

                        local_test.swap_once();
                        double new_score = local_test.get_normality_statistic();
                        if (new_score < score) {
                            if (new_score < best_thread.fitness()) {
                                best_thread = solution(local_test.current_subset(),
                                                       local_test.get_normality_statistic());
                            }
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
                    return best_thread;
                }));
            }

            solution the_best;
            for (auto& future: thread_solutions) {
                the_best = std::min(the_best, future.get());
            }
            best.push_back(the_best.elems());
            best_scores.push_back(the_best.fitness());

            test.set_solution(the_best.elems());

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

