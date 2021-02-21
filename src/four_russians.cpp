#include "include/mvn_test.h"

mvn::four_russians::four_russians(const mvn::Matrix& pairwise_mahalanobis, size_t batch_size) : current(pairwise_mahalanobis.size(),
                                                                                                        batch_size) {
    n = pairwise_mahalanobis.size();
    for (int i = 0; i < n; i++) {
        batches.emplace_back();
        // round up division
        unsigned number_of_states = (n + batch_size - 1) / batch_size;
        for (int j = 0; j < number_of_states; j++) {
            batches.back().emplace_back(pairwise_mahalanobis.row(i), j * batch_size, std::min((j + 1) * batch_size, n));
        }
    }
}

void mvn::four_russians::set(unsigned int i) {
    current.set(i);
}

void mvn::four_russians::unset(unsigned int i) {
    current.unset(i);
}

double mvn::four_russians::score(unsigned int row) const {
    auto states = current.get_states();
    double score = 0.0;
    for (unsigned i = 0; i < batches.size(); i++) {
        score += batches[row][i].get_batch_score(states[i]);
    }
    return score;
}

bool mvn::four_russians::is_effective() const {
    return current.effective_size() > batches.size() * 2;
}

mvn::batch::batch(const mvn::Vector& v, unsigned int start, unsigned int stop) {
    scores.push_back(0.0);
    unsigned n = stop - start;
    for (unsigned i = 0; i < (1 << n); i++) {
        double score = 0.0;
        for (unsigned j = 0; j < n; j++) {
            if (i & (1 << j)) {
                score += v(start + j);
            }
        }
        scores.push_back(score);
    }
}

double mvn::batch::get_batch_score(unsigned int state) const {
    return scores.at(state);
}

mvn::state::state(size_t size, size_t chunk_size) : chunk_size(chunk_size), subset_size(0) {
    bits.resize(size / chunk_size, 0);
}
