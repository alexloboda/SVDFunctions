#include "include/RandomSampler.h"
#include <cmath>
#include <cassert>
#include <limits>

namespace mvn {

RandomSampler::RandomSampler() : runif(0.0, 1.0), wheel(0), size(0) {}

RandomSampler::RandomSampler(const std::vector<double>& logscale, long seed)
    : runif(0.0, 1.0), wheel(seed), original(logscale), size(logscale.size()) {
    if (original.size() < 2) {
        throw std::invalid_argument("Too little segment tree");
    }

    active_tree.resize(2 * original.size() - 1);
    segment_tree.resize(2 * original.size() - 1);
    for (size_t i = 0; i < original.size(); i++) {
        auto k = segment_tree.size() - i - 1;
        active_tree[k] = 1;
        segment_tree[k] = original[i];
    }
    for (int i = segment_tree.size() - size - 1; i >= 0; i--) {
        update_inner_node(i);
    }
}

RandomSampler::RandomSampler(const RandomSampler& other)
    : runif(0.0, 1.0), wheel(other.wheel), original(other.original),
      segment_tree(other.segment_tree), active_tree(other.active_tree), size(other.size) {}

RandomSampler& RandomSampler::operator=(RandomSampler&& other) {
    runif = other.runif;
    wheel = other.wheel;
    original = std::move(other.original);
    segment_tree = std::move(other.segment_tree);
    active_tree = std::move(other.active_tree);
    size = other.size;
    return *this;
}

void RandomSampler::disable(size_t n) {
    if (!is_active(n)) {
        throw std::logic_error("Disabling non-active element");
    }
    auto pos = el_pos(n);
    segment_tree[pos] = -std::numeric_limits<double>::infinity();
    active_tree[pos] = false;
    update(pos);
}

void RandomSampler::enable(size_t n) {
    if (is_active(n)) {
        throw std::logic_error("Enabling active element.");
    }
    auto pos = el_pos(n);
    segment_tree[pos] = original[n];
    active_tree.at(pos) = true;
    update(pos);
}

size_t RandomSampler::sample() {
    if (size < 2) {
        throw std::logic_error("Too little segment tree");
    }
    size_t node = 0;
    while (!is_leaf(node)) {
        auto chld = children(node);
        if (active_tree.at(chld.first) == 0) {
            if (active_tree.at(chld.second) == 0) {
                throw std::logic_error("No active elements in subtree");
            }
            node = chld.second;
        } else if (active_tree.at(chld.second) == 0) {
            node = chld.first;
        } else {
            double l = segment_tree.at(chld.first);
            double r = segment_tree.at(chld.second);
            double maxL = std::max(l, r);
            l = std::exp(l - maxL);
            r = std::exp(r - maxL);
            auto sum = l + r;
            l /= sum;
            if (runif(wheel) < l) {
                node = chld.first;
            } else {
                node = chld.second;
            }
        }
    }

    int aug_nodes = segment_tree.size() - original.size();
    if (!active_tree.at(node)) {
        throw std::logic_error("Sampled element is not active.");
    }
    int element = node - aug_nodes;
    assert(el_pos(element) == node);
    return element;
}

size_t RandomSampler::n_active() const {
    return active_tree.at(0);
}

std::pair<size_t, size_t> RandomSampler::children(size_t node) const {
    if (node > segment_tree.size() - original.size() - 1) {
        throw std::invalid_argument("It's a leaf");
    }
    return {2 * node + 1, 2 * node + 2};
}

bool RandomSampler::is_root(size_t node) {
    return node == 0;
}

bool RandomSampler::is_leaf(size_t node) const {
    return node >= segment_tree.size() - original.size();
}

bool RandomSampler::is_active(size_t node) const {
    return active_tree.at(el_pos(node));
}

size_t RandomSampler::el_pos(size_t el) const {
    size_t pos = segment_tree.size() - original.size() + el;
    if (pos >= segment_tree.size() || pos < 0) {
        throw std::out_of_range("Out of range");
    }
    return pos;
}

double RandomSampler::sum_log(double l, double r) {
    if (std::isinf(l)) {
        return r;
    }
    if (std::isinf(r)) {
        return l;
    }
    double maxL = std::max(l, r);
    l = std::exp(l - maxL);
    r = std::exp(r - maxL);
    return std::log(l + r) + maxL;
}

size_t RandomSampler::parent(size_t node) {
    return (node - 1) / 2;
}

void RandomSampler::update_inner_node(size_t node) {
    auto chs = children(node);
    segment_tree[node] = sum_log(segment_tree[chs.first], segment_tree[chs.second]);
    assert(!std::isnan(segment_tree[node]));
    active_tree[node] = active_tree.at(chs.first) + active_tree.at(chs.second);
}

void RandomSampler::update(size_t node) {
    if (node < segment_tree.size() - original.size()) {
        throw std::invalid_argument("That's not a leaf");
    }
    node = parent(node);
    while (true) {
        update_inner_node(node);
        if (is_root(node)) {
            break;
        }
        node = parent(node);
    }
}

} // namespace mvn
