#ifndef SRC_VCF_PREDICTING_HANDLER_H
#define SRC_VCF_PREDICTING_HANDLER_H

#include "vcf_handlers.h"

namespace vcf {
    class PredictingHandler : public VariantsHandler {
        const int window_size_kb;
        const int window_size;

        GenotypeMatrixHandler& gh;
        Chromosome curr_chr;
        std::unique_ptr<Variant> awaiting;
        std::unordered_map<int, std::set<Range>> ranges;
    public:
        explicit PredictingHandler(const std::vector<std::string>& samples, GenotypeMatrixHandler& gh,
                                   int window_size_kb, int window_size);
        void processVariant(const Variant& variant, const std::vector<Allele>& alleles) override;
        bool isOfInterest(const Variant& position) override;
        void cleanup();
    };
}

#endif //SRC_VCF_PREDICTING_HANDLER_H
