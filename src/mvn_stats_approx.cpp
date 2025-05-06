#include "include/mvn_stats_approx.h"
#include "include/kronecker.h"
#include "include/mahalanobis_distances.h"
#include "include/mvn_clst.h"

namespace mvn {

mvn_stats_approx::mvn_stats_approx(const std::string& _filename, int n_clusters)
    : mvn_stats(n_clusters), filename(_filename) {}

void mvn_stats_approx::init_pairwise(mahalanobis_distances distances, const Clustering& clst, double beta) {
    double k_pw = -(beta * beta) / 2.0;

    auto n = clst.size();
    matching::kronecker_calculator calc(filename);
    mahalanobis_pairwise = calc.calculate(distances.get_sigma(), k_pw);

    for (size_t cl = 0; cl < n; cl++) {
        for (size_t pair_cl = 0; pair_cl < n; pair_cl++) {
            auto size_cl = clst.elements(cl).size();
            auto size_pair_cl = clst.elements(pair_cl).size();
            double bias = (size_cl * size_pair_cl);
            mahalanobis_pairwise[cl][pair_cl] += bias;
        }
    }

    for (size_t cl = 0; cl < n; cl++) {
        mahalanobis_pairwise[cl][cl] *= 0.5;
    }
}

} // namespace mvn
