#include <vector>

namespace clustering {

double euclidean_distance(const std::vector<double>& x, const std::vector<double>& y);

class Point {
public:
    int idx;
    int current_assignment;
    double dist;

    Point(int _idx, int _current_assignment, double _dist);
    bool operator<(const Point& other) const;
};

class SameSizeKMeans {
    int k, max_size, max_iter;
    double tol;
    std::vector<std::vector<double>> centroids;
    std::vector<int> labels;
    std::vector<int> cluster_size;

    std::vector<std::vector<double>> k_means_plus_plus(const std::vector<std::vector<double>>& X, int k);
    void k_means(const std::vector<std::vector<double>>& X, int k); 

    void initialize(const std::vector<std::vector<double>>& X);    
    void update_centroids(const std::vector<std::vector<double>>& X);
    double inertia(const std::vector<std::vector<double>>& X) const;
public:
    SameSizeKMeans(int k, int max_iter = 10000, double tol = 1e-4);
    void fit(const std::vector<std::vector<double>>& X);
    std::vector<int> get_labels();
};

}