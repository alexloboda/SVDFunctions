#include "include/vcf_predicting_handler.h"

namespace vcf {
    void insert(Range r, std::set<Range>& ranges) {
        auto pos = ranges.lower_bound(r);
        if (pos != ranges.end() && pos->includes(r.end() + (-1))) {
            int from = std::min(r.begin().position(), pos->begin().position());
            int to = std::max(r.end().position(), pos->end().position());
            r = Range(r.begin().chromosome(), from, to);
            ranges.erase(*pos);
        }
        ranges.insert(r);
    }

    PredictingHandler::PredictingHandler(const std::vector<std::string>& samples, GenotypeMatrixHandler& gh,
                                         int window_size_kb, int window_size)
                                         :VariantsHandler(samples), window_size_kb(window_size_kb),
                                         window_size(window_size), gh(gh), curr_chr(-1), awaiting(nullptr) {
        auto variants = gh.desired_variants();
        int halfws = window_size_kb / 2;
        for (const Variant& v : variants) {
            Chromosome chr = v.position().chromosome();
            int pos = v.position().position();
            Range haplotype{chr, pos - halfws, pos + halfws};
            insert(haplotype, ranges[chr.num()]);
        }
    }

    bool PredictingHandler::isOfInterest(const Variant& variant) {
        Position pos = variant.position();
        // Workaround, c++11
        auto& set = ranges[pos.chromosome().num()];
        auto it = set.lower_bound(Range(pos.chromosome(), -1, pos.position()));
        return it != set.end() && it->includes(pos);
    }
}
