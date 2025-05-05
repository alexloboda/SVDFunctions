#ifndef MVN_STATS_APPROX_H
#define MVN_STATS_APPROX_H

#include "mvn_stats.h"

namespace mvn {

class mvn_stats_approx : public mvn_stats {
public:
    mvn_stats_approx(const std::string& filename, int n_clusters);
    void init_pairwise(mahalanobis_distances distances, const Clustering& clst, double beta) override;

private:
    std::string filename;
};

} // namespace mvn

#endif // MVN_STATS_APPROX_H
