#include "include/mvn_stats_interpoint.h"
#include "include/mahalanobis_distances.h"
#include "include/mvn_clst.h"

namespace mvn {

mvn_stats_interpoint::mvn_stats_interpoint(int n_clusters) : mvn_stats(n_clusters) {}

void mvn_stats_interpoint::init_pairwise(mahalanobis_distances distances, const Clustering& clst, double beta) {
    distances.calculate_interpoint();

    auto n = clst.size();
    double k_pw = -(beta * beta) / 2.0;

    auto calculate_contribution = [&](size_t cl, size_t pair_cl) {
        double contribution = 0.0;
        for (int el : clst.elements(cl)) {
            for (int pair_el : clst.elements(pair_cl)) {
                double distance = distances.interpoint_distance(el, pair_el);
                double weight = (cl == pair_cl) ? 0.5 : 1.0;
                double additive = std::exp(k_pw * distance);
                contribution += weight * additive;
            }
        }
        return contribution;
    };

    for (size_t cl = 0; cl < n; cl++) {
        for (size_t pair_cl = 0; pair_cl < n; pair_cl++) {
            mahalanobis_pairwise[cl][pair_cl] = calculate_contribution(cl, pair_cl);
        }
    }
}

} // namespace mvn
