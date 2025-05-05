#ifndef RANDOM_SAMPLER_H
#define RANDOM_SAMPLER_H

#include <vector>
#include <random>
#include <stdexcept>

namespace mvn {

class RandomSampler {
    std::uniform_real_distribution<double> runif;
    mutable std::mt19937 wheel;
    std::vector<double> original;
    std::vector<double> segment_tree;
    std::vector<size_t> active_tree;

    size_t size;

public:
    RandomSampler();
    RandomSampler(const std::vector<double>& logscale, long seed);
    RandomSampler(RandomSampler&&) = default;
    RandomSampler(const RandomSampler& other);
    RandomSampler& operator=(RandomSampler&&);

    void disable(size_t n);
    void enable(size_t n);
    size_t sample();
    size_t n_active() const;

private:
    std::pair<size_t, size_t> children(size_t node) const;
    static bool is_root(size_t node);
    bool is_leaf(size_t node) const;
    bool is_active(size_t node) const;
    size_t el_pos(size_t el) const;
    static size_t parent(size_t node);

    static double sum_log(double l, double r);

    void update_inner_node(size_t node);
    void update(size_t node);
};

} // namespace mvn

#endif // RANDOM_SAMPLER_H
