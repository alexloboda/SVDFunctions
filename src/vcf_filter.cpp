#include <algorithm>
#include "include/vcf_filter.h"

namespace {
    using std::vector;
    using std::string;
}

namespace vcf {
    VCFFilter::VCFFilter(int DP, int GQ)
        :samples_set(false), DP(DP), GQ(GQ) {}

    void VCFFilter::add_bad_variants(vector<Position>& positions) {
        bad_positions.insert(positions.begin(), positions.end());
    }

    void VCFFilter::add_samples(vector<string>& samples) {
        samples_set = true;
        available_samples.insert(samples.begin(), samples.end());
    }

    bool VCFFilter::apply(const Position& p) const {
        return bad_positions.find(p) == bad_positions.end();
    }

    bool VCFFilter::apply(const string& sample) const {
        return !samples_set || available_samples.find(sample) != available_samples.end();
    }

    bool VCFFilter::apply(int dp, int gq) const {
        return dp >= DP && gq >= GQ;
    }
}
