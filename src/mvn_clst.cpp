#include "include/mvn_clst.h"

namespace mvn {

Clustering::Clustering(const std::vector<int>& clustering) {
    if (clustering.empty()) {
        throw std::invalid_argument("Clustering must not be empty");
    }

    int n_clusters = *std::max_element(clustering.begin(), clustering.end()) + 1;
    cluster_sizes.resize(n_clusters);
    clusters.resize(n_clusters);
    for (size_t i = 0; i < clustering.size(); i++) {
        int cl = clustering[i];
        ++cluster_sizes[cl];
        clusters[cl].push_back(i);
    }
}

const std::vector<int>& Clustering::elements(size_t i) const {
    return clusters.at(i);
}

size_t Clustering::size() const {
    return clusters.size();
}

size_t Clustering::cluster_size(size_t i) const {
    return clusters.at(i).size();
}

} // namespace mvn
