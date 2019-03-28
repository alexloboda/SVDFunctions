#ifndef SRC_VCF_STATS_H
#define SRC_VCF_STATS_H

#include <vector>
#include <limits>
#include <string>

namespace vcf {
    enum class Stat {
        OVERALL, NON_PASS, BANNED, WARNING, MISSING_RATE, GT_MISS, DP_GQ, ALLELE_BALANCE
    };

    std::vector<Stat> statsList();
    std::string to_string(Stat stat);

    class VCFFilterStats {
        std::vector<unsigned long> stats;
    public:
        VCFFilterStats();
        void add(Stat, unsigned long);
        unsigned long value(Stat);
    };
}
#endif //SRC_VCF_STATS_H
