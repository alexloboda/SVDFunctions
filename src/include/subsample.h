#ifndef SRC_SUBSAMPLE_H
#define SRC_SUBSAMPLE_H

#include <vector>
#include "mvn_test.h"
namespace mvn {

class subsample {
    std::shared_ptr<mvn_test> test;

    std::vector<std::vector<size_t>> best;
    Clustering clst;

    std::mt19937 wheel;
public:
    subsample();
    subsample(std::shared_ptr<const Matrix> X, const Clustering& clst, const Vector& mean, const Matrix& cov);
    subsample(subsample&&) = default;
    subsample& operator=(subsample&& other) = default;

    void run(size_t iterations, size_t restarts, double t0, double c, size_t threads, size_t start, size_t size_ub, size_t step);

    std::vector<size_t> get_solution(size_t size) const;
    size_t solutions() const;
};

}
#endif //SRC_SUBSAMPLE_H
