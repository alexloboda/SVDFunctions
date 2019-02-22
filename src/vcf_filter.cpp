#include <algorithm>
#include "vcf_filter.h"

namespace {
    using std::vector;
    using std::string;
}

namespace vcf {
    VCFFilter::VCFFilter(int DP, int GQ)
        :variants_set(false), samples_set(false), DP(DP), GQ(GQ) {}

    void VCFFilter::set_available_variants(vector<Variant>& variants) {
        available_variants.clear();
        available_variants.insert(variants.begin(), variants.end());
        variants_set = true;
    }

    void VCFFilter::add_bad_variants(vector<Position>& positions) {
        bad_variants.insert(positions.begin(), positions.end());
    }

    void VCFFilter::add_samples(vector<string>& samples) {
        samples_set = true;
        available_samples.insert(samples.begin(), samples.end());
    }

    bool VCFFilter::apply(const vcf::Variant& v) const {
        return (!variants_set || available_variants.find(v) != available_variants.end());
    }

    bool VCFFilter::apply(const Position& p) const {
        return bad_variants.find(p) == bad_variants.end();
    }

    bool VCFFilter::apply(const string& sample) const {
        return !samples_set || available_samples.find(sample) != available_samples.end();
    }

    bool VCFFilter::apply(int dp, int gq) const {
        return dp > DP && gq > GQ;
    }
}
