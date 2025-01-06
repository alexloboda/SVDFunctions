#include "include/sskm.h"

#include <cmath>
#include <queue>
#include <random>
#include <algorithm>
#include <iostream>
#include <assert.h>

#include <Rcpp.h>

namespace {
    class Preference {
    public:
        int idx;
        int cluster;
        double value;

        Preference(int _idx, int _cluster, double _value) {
            idx = _idx;
            cluster = _cluster;
            value = _value;
        }

        bool operator<(const Preference& other) const {
            return value < other.value;
        }
    };
}

namespace clustering {

 
Point::Point(int _idx, int _current_assignment, double _dist) {
    idx = _idx;
    current_assignment = _current_assignment;
    dist = _dist;
}

bool Point::operator<(const Point& other) const {
    return dist < other.dist;
}

double euclidean_distance(const std::vector<double>& x, const std::vector<double>& y) {
    double sum = 0;
    for (int i = 0; i < x.size(); i++) {
        sum += (x[i] - y[i]) * (x[i] - y[i]);
    }
    return std::sqrt(sum);

}

SameSizeKMeans::SameSizeKMeans(int _k, int _max_iter, double _tol) {
    k = _k;
    max_iter = _max_iter;
    tol = _tol;
}

std::vector<std::vector<double>> SameSizeKMeans::k_means_plus_plus(const std::vector<std::vector<double>>& X, int k) {
    std::vector<std::vector<double>> centroids;
    std::vector<double> dists(X.size(), 1e9);
    std::uniform_int_distribution<int> dist(0, X.size() - 1);
    std::mt19937 gen(42);
    int idx = dist(gen);
    centroids.push_back(X[idx]);
    for (int i = 1; i < k; i++) {
        for (int j = 0; j < X.size(); j++) {
            double d = euclidean_distance(X[j], centroids[i - 1]);
            d = d * d;
            if (d < dists[j]) {
                dists[j] = d;
            }
        }
        std::discrete_distribution<int> dist(dists.begin(), dists.end());
        idx = dist(gen);
        centroids.push_back(X[idx]);
    }
    return centroids;
}

std::vector<int> SameSizeKMeans::get_labels() {
    return labels;
}

void SameSizeKMeans::initialize(const std::vector<std::vector<double>>& X) {
    int n = X.size();
    int assigned = 0;
    std::vector<int> centroids_left;
    std::vector<int> points_to_assign;

    for (int i = 0; i < n; i++) {
        points_to_assign.push_back(i);
    }

    for (int i = 0; i < centroids.size(); i++) {
        centroids_left.push_back(i);
    }

    while (assigned < n) {
        std::vector<Preference> preferences;

        for (int i: points_to_assign) {
            std::vector<double> dists(centroids_left.size());
            for (int j = 0; j < centroids_left.size(); j++) {
                dists[j] = euclidean_distance(X[i], centroids[centroids_left[j]]);
            }
            auto min_idx = std::min_element(dists.begin(), dists.end()) - dists.begin();
            auto max_idx = std::max_element(dists.begin(), dists.end()) - dists.begin();

            double score = dists[min_idx] - dists[max_idx];
            Preference p(i, centroids_left.at(min_idx), score);
            preferences.push_back(p);
        }

        std::sort(preferences.begin(), preferences.end());

        for (Preference p: preferences) {
            int idx = p.idx;
            int candidate = p.cluster;

            if (cluster_size[candidate] < max_size) {
                labels[idx] = candidate;
                cluster_size[candidate]++;
                assigned++;
            } else {
                centroids_left.erase(std::remove(centroids_left.begin(), centroids_left.end(), candidate), centroids_left.end());
                auto points_to_assign_new = std::vector<int>();
                for (int i: points_to_assign) {
                    if (labels[i] == -1) {
                        points_to_assign_new.push_back(i);
                    }
                }
                points_to_assign = points_to_assign_new;
                break;
            }

        }
    }

    for (int i = 0; i < n; i++) {
        assert(labels[i] != -1);
    }

    return;
}

namespace {
    void update_heap(int point, int new_cluster, 
                    const std::vector<std::vector<double>>& X,
                    const std::vector<int> labels, 
                    const std::vector<std::vector<double>>& centroids, 
                    std::vector<std::vector<std::priority_queue<Point>>>& heap, 
                    double tol) {
        int k = centroids.size();
        for (int i = 0; i < k; i++) {
            if (i == new_cluster) {
                continue;
            }
            while (!heap[new_cluster][i].empty()) {
                Point p = heap[new_cluster][i].top();
                if (p.current_assignment != labels[p.idx]) {
                    heap[new_cluster][i].pop();
                } else {
                    break;
                }
            }

            double dist_cluster = euclidean_distance(X[point], centroids[i]);
            double current_dist = euclidean_distance(X[point], centroids[new_cluster]);
            double gain = current_dist - dist_cluster;
            Point p(point, new_cluster, gain);
            heap[new_cluster][i].push(p);
        }
    }
}

double SameSizeKMeans::inertia(const std::vector<std::vector<double>>& X) const {
    double sum = 0;
    for (int i = 0; i < X.size(); i++) {
        sum += euclidean_distance(X[i], centroids[labels[i]]);
    }
    return sum;
}

void SameSizeKMeans::update_centroids(const std::vector<std::vector<double>>& X) {
    int n = X.size();
    int k = centroids.size();
    std::vector<std::vector<double>> new_centroids(k, std::vector<double>(X[0].size(), 0));
    std::vector<int> cluster_size(k, 0);

    for (int i = 0; i < n; i++) {
        int cluster = labels[i];
        for (int j = 0; j < X[i].size(); j++) {
            new_centroids[cluster][j] += X[i][j];
        }
        cluster_size[cluster]++;
    }

    for (int i = 0; i < k; i++) {
        for (int j = 0; j < X[0].size(); j++) {
            new_centroids[i][j] /= cluster_size[i];
        }
    }

    centroids = new_centroids;
}

void SameSizeKMeans::k_means(const std::vector<std::vector<double>>& X, int k) {
    centroids = k_means_plus_plus(X, k);
    labels.resize(X.size(), -1);
    int n = X.size();
    cluster_size.resize(k, 0);

    initialize(X);
    update_centroids(X);
    
    for (int iter = 0; iter < max_iter; iter++) {
        double prev_inertia = inertia(X);
        std::vector<std::vector<std::priority_queue<Point>>> closest_transfers(k);

        closest_transfers.resize(k);
        for (int i = 0; i < k; i++) {
            closest_transfers[i].resize(k);
        }

        std::priority_queue<Point> points_order;
        std::vector<double> best_dist(n, 0);
        for (int i = 0; i < n; i++) {
            int cluster = labels[i];
            std::vector<double> dists(k);
            
            int min_idx = -1;
            for (int j = 0; j < k; j++) {
                dists[j] = euclidean_distance(X[i], centroids[j]);
                if (j != cluster) {
                    if (min_idx == -1 || dists[j] < dists[min_idx]) {
                        min_idx = j;
                    }
                }
            }

            best_dist[i] = dists[min_idx];
            double priority = dists[labels[i]] - dists[min_idx];
            Point p(i, cluster, priority);   
            points_order.push(p);
        }

        double total_gain = 0;

        while (!points_order.empty()) {
            Point p = points_order.top();
            points_order.pop();
            int idx = p.idx;
            int cluster = labels[idx];
            int priority = p.dist;

            std::vector<double> dists(k);
            for (int j = 0; j < k; j++) {
                dists[j] = euclidean_distance(X[idx], centroids[j]);
            }

            auto min_idx = std::min_element(dists.begin(), dists.end()) - dists.begin();
            auto min_dist = dists[min_idx];

            // Check if inertia can be improveed
            if (std::fabs(min_dist - dists[cluster]) < tol) {
                continue;
            }

            int best_cluster = -1;
            double max_gain = std::numeric_limits<double>::lowest();
            bool moved = false;

            for (int j = 0; j < k; j++) {
                if (j == cluster) {
                    continue;
                }   
    
                double gain = dists[cluster] - dists[j];
                while (!closest_transfers[j][cluster].empty()) {
                    Point closest = closest_transfers[j][cluster].top();
                    if (closest.current_assignment != labels[closest.idx]) {
                        closest_transfers[j][cluster].pop();
                        continue;
                    }
                    
                    bool this_moved = false;
                    if (closest.dist < 0 && cluster_size[j] < max_size) {
                        this_moved = true;  
                    } else {
                        gain += closest.dist;
                    }

                    assert(gain <= prev_inertia);
                    if (gain > max_gain && gain > tol) {
                        max_gain = gain;
                        best_cluster = j;
                        moved = this_moved;
                    }
                    break;
                }
            }

            if (best_cluster == -1) {
                // add current point to transer list
                for (int i = 0; i < k; i++) {
                    if (i != cluster) {
                        double gain = dists[cluster] - dists[i];
                        Point p(idx, cluster, gain);
                        closest_transfers[cluster][i].push(p);
                    }
                }
                continue;
            } 

            total_gain += max_gain;

            if (moved) {
               labels[idx] = best_cluster;
               cluster_size[best_cluster]++;
               cluster_size[cluster]--; 
               continue;
            }

            Point candidate = closest_transfers[best_cluster][cluster].top();
            closest_transfers[best_cluster][cluster].pop();
            labels[idx] = best_cluster;
            labels[candidate.idx] = cluster;

            update_heap(candidate.idx, cluster, X, labels, 
                        centroids, closest_transfers, tol);
            update_heap(idx, best_cluster, X, labels, 
                        centroids, closest_transfers, tol);
        }

        double current_inertia = inertia(X);
        assert (fabs(prev_inertia - current_inertia - total_gain) <= tol);
        update_centroids(X);
        std::cout << current_inertia << std::endl;
        assert (current_inertia - prev_inertia < tol);

        if (prev_inertia - current_inertia < tol) {
            break;
        }
        prev_inertia = current_inertia;
    }
}
    
void SameSizeKMeans::fit(const std::vector<std::vector<double>>& X) {
    cluster_size.clear();
    labels.clear();
    max_size = X.size() / k;
    if (X.size() % k != 0) {
        max_size++;
    }       
    k_means(X, k);

}

}