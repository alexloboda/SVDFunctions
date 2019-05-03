#include "include/vcf_stats.h"
#include <stdexcept>

namespace vcf {

    VCFFilterStats::VCFFilterStats() {
        stats.resize(statsList().size());
    }

    void VCFFilterStats::add(Stat stat, unsigned long num) {
        stats.at((unsigned)stat) += num;
    }

    unsigned long VCFFilterStats::value(Stat stat) {
        return stats.at((unsigned)stat);
    }

    std::vector<Stat> statsList() {
        return {Stat::OVERALL, Stat::NON_PASS, Stat::BANNED, Stat::WARNING, Stat::MISSING_RATE,
                Stat::GT_MISS, Stat::DP_GQ, Stat::ALLELE_BALANCE};
    }

    std::string to_string(Stat stat) {
        switch (stat) {
            case Stat::OVERALL:
                return "OVERALL";
            case Stat::NON_PASS:
                return "NON_PASS";
            case Stat::BANNED:
                return "BANNED";
            case Stat::WARNING:
                return "WARNING";
            case Stat::MISSING_RATE:
                return "MISSING_RATE";
            case Stat::GT_MISS:
                return "GT_MISS";
            case Stat::DP_GQ:
                return "GT_DP_GQ_FILTERS";
            case Stat::ALLELE_BALANCE:
                return "Allele balance";
            default:
                throw std::logic_error("Unreachable statement");
        }
    }

}
