#ifndef SRC_VCF_FILTER_H
#define SRC_VCF_FILTER_H

#include "vcf_primitives.h"

namespace vcf {

    class VCFFilter {
        int DP;
        int GQ;

        std::unordered_set<Position> bad_positions;
        std::unordered_set<std::string> available_samples;

        bool samples_set;
    public:
        VCFFilter(int DP, int GQ);
        void add_bad_variants(std::vector<Position>& positions);
        void add_samples(std::vector<std::string>& samples);
        bool apply(const Position& v) const;
        bool apply(int dp, int gq) const;
        bool apply(const std::string& sample) const;
    };

}

#endif //SRC_VCF_FILTER_H
