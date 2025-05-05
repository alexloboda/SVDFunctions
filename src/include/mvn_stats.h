#ifndef MVN_STATS_H
#define MVN_STATS_H

#include <vector>
#include <memory>
#include "mahalanobis_distances.h"
#include "mvn_clst.h"

namespace mvn {

class mvn_stats {
public:
    explicit mvn_stats(int n_clusters);
    virtual void init_pairwise(mahalanobis_distances distances, const Clustering& clst, double beta) = 0;
    void init(mahalanobis_distances distances, const Clustering& clst, double beta); 

    double pairwise_stat(size_t i, size_t j) const;
    double centered_stat(size_t i) const;
    double sum_pairwise(size_t point, const std::vector<size_t>& ss) const;

protected:
    std::vector<std::vector<double>> mahalanobis_pairwise;
    std::vector<double> mahalanobis_centered;
};

} // namespace mvn

#endif // MVN_STATS_H
