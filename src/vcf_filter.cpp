#include <algorithm>
#include "vcf_filter.h"

namespace {
    using std::vector;
    using std::string;
}

namespace vcf {
    VCFFilter::VCFFilter(int DP, int GQ)
        :variants_set(false), DP(DP), GQ(GQ) {}

    void VCFFilter::set_available_variants(vector<Variant>& variants) {
        available_variants.clear();
        available_variants.insert(variants.begin(), variants.end());
        variants_set = true;
    }

    void VCFFilter::add_bad_variants(vector<Variant>& variants) {
        bad_variants.insert(variants.begin(), variants.end());
    }

    void VCFFilter::add_samples(vector<string> samples) {
        available_samples.insert(samples.begin(), samples.end());
    }

    bool VCFFilter::apply(const Variant& v) const {
        return bad_variants.find(v) == bad_variants.end() &&
                (!variants_set || available_variants.find(v) != available_variants.end());
    }

    bool VCFFilter::apply(const string& sample) const {
        return available_samples.find(sample) != available_samples.end();
    }

    bool VCFFilter::apply(const Allele& allele) const {
        return allele.DP() > DP && allele.GQ() > GQ;
    }
}
