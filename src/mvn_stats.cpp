#include "include/mvn_stats.h"

namespace mvn {

mvn_stats::mvn_stats(int n_clusters) :mahalanobis_centered(n_clusters), mahalanobis_pairwise(n_clusters) {
    for (size_t cl = 0; cl < n_clusters; cl++) {
        mahalanobis_pairwise[cl].resize(n_clusters);
    }
}

void mvn_stats::init(mahalanobis_distances distances, const Clustering& clst, double beta) {
    size_t n = clst.size();
    double k_center = -(beta * beta / (2 * (1 + beta * beta)));

    for (size_t cl = 0; cl < n; cl++) {
        for (int el : clst.elements(cl)) {
            mahalanobis_centered[cl] += std::exp(k_center * distances.distance(el));
        }
    }
}

double mvn_stats::pairwise_stat(size_t i, size_t j) const {
    return mahalanobis_pairwise[i][j];
}

double mvn_stats::centered_stat(size_t i) const {
    return mahalanobis_centered[i];
}

double mvn_stats::sum_pairwise(size_t point, const std::vector<size_t>& ss) const {
    double ret = 0.0;
    for (auto s: ss) {
        ret += mahalanobis_pairwise[point][s];
    }
    return ret;
}

} // namespace mvn