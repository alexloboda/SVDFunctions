#ifndef SRC_SUBSAMPLE_H
#define SRC_SUBSAMPLE_H

#include <vector>
#include "mvn_test.h"
namespace mvn {
class subsample {
    std::shared_ptr<mvn_test> test;

    std::vector<std::vector<size_t>> best;
    std::vector<double> best_scores;
    Clustering clst;

    std::mt19937 wheel;
public:
    subsample();
    subsample(std::shared_ptr<const Matrix> X, const Clustering& clst, const Vector& mean, const Matrix& cov);
    subsample(std::shared_ptr<const Matrix> X, const Clustering& clst);
    subsample(subsample&&) = default;
    subsample& operator=(subsample&& other) = default;

    void run(size_t iterations, size_t restarts, double t0, int threads, int start, int size_ub, int step);

    std::vector<size_t> get_solution(size_t size) const;
    size_t solutions() const;
};

}
#endif //SRC_SUBSAMPLE_H
