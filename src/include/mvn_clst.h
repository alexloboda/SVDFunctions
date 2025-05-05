#ifndef CLUSTERING_H
#define CLUSTERING_H

#include <vector>
#include <stdexcept>
#include <algorithm>

namespace mvn {

class Clustering {
    std::vector<int> cluster_sizes;
    std::vector<std::vector<int>> clusters;

public:
    explicit Clustering(const std::vector<int>& clustering);
    Clustering() = default;

    const std::vector<int>& elements(size_t i) const;
    size_t size() const;
    size_t cluster_size(size_t i) const;
};

} // namespace mvn

#endif // CLUSTERING_H
