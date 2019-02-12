#ifndef SRC_VCF_FILTER_H
#define SRC_VCF_FILTER_H

#include "vcf_primitives.h"

namespace vcf {

    class VCFFilter {
        int DP;
        int GQ;

        std::unordered_set<Variant> available_variants;
        std::unordered_set<Variant> bad_variants;
        std::unordered_set<std::string> available_samples;
        std::vector<Range> ranges;

        bool ranges_set;
        bool variants_set;
    public:
        VCFFilter(int DP, int GQ);
        void set_available_variants(std::vector<Variant>& variants);
        void add_bad_variants(std::vector<Variant>& variants);
        void set_ranges(std::vector<Range>& ranges);
        void add_samples(std::vector<std::string> samples);
        bool apply(const Variant& v) const;
        bool apply(const std::string& sample) const;
        bool apply(const Allele& allele) const;
    };

}

#endif //SRC_VCF_FILTER_H
