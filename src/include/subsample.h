#ifndef SRC_SUBSAMPLE_H
#define SRC_SUBSAMPLE_H

#include <vector>
#include "mvn_test.h"
namespace mvn {

class subsample {
    std::shared_ptr<mvn_test> test;

    std::vector<std::vector<size_t>> best;
    std::vector<double> best_stat;

    // Sole source of randomness for the whole search: every worker's stream is seeded
    // from here, on the dispatching thread. Same seed, same build, same machine ->
    // same answer; see the note on clone() in mvn_test.h.
    std::mt19937 wheel;
public:
    subsample();
    subsample(std::shared_ptr<const Matrix> X, const Clustering& clst, const Vector& mean, const Matrix& cov,
              const PrecomputeConfig& config, std::mt19937::result_type seed);
    subsample(subsample&&) = default;
    subsample& operator=(subsample&& other) = default;

    void run(size_t iterations, size_t restarts, double t0, double c, size_t threads, size_t start, size_t size_ub,
             size_t step);

    std::vector<size_t> get_solution(size_t size) const;
    size_t solutions() const;
    double statistic(size_t k);
};

}
#endif //SRC_SUBSAMPLE_H
