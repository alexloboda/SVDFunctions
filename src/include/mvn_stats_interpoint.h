#ifndef MVN_STATS_INTERPOINT_H
#define MVN_STATS_INTERPOINT_H

#include "mvn_stats.h"

namespace mvn {

class mvn_stats_interpoint : public mvn_stats {
public:
    explicit mvn_stats_interpoint(int n_clusters);
    void init_pairwise(mahalanobis_distances distances, const Clustering& clst, double beta) override;
};

} // namespace mvn

#endif // MVN_STATS_INTERPOINT_H
