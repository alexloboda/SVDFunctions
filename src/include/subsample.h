#ifndef SRC_SUBSAMPLE_H
#define SRC_SUBSAMPLE_H

#include <vector>
#include "mvn_test.h"
namespace mvn {
class subsample {
    mvn_test test;

    std::vector<std::vector<size_t>> best;
    std::vector<double> best_scores;

    std::mt19937 wheel;
    std::uniform_real_distribution<double> random_unif;
public:
    subsample();
    subsample(const Matrix& X, const Vector& mean, const Matrix& cov);
    subsample(const subsample&) = default;
    subsample& operator=(const subsample& other) = default;

    std::vector<size_t> run(size_t iterations, size_t restarts, double t0);

    std::vector<size_t> get_solution(size_t size) const;
    double get_statistic(size_t size) const;
    size_t min_size() const;
    size_t max_size() const;
};

}
#endif //SRC_SUBSAMPLE_H
