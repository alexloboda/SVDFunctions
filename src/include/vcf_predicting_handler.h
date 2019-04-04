#ifndef SRC_VCF_PREDICTING_HANDLER_H
#define SRC_VCF_PREDICTING_HANDLER_H

#include "vcf_handlers.h"

namespace vcf {
    class PredictingHandler : public VariantsHandler {
        const unsigned long window_size_kb;
        const unsigned long window_size;

        GenotypeMatrixHandler& gh;
        Chromosome curr_chr;
        std::unique_ptr<Variant> awaiting;

    public:
        explicit PredictingHandler(const std::vector<std::string>& samples, GenotypeMatrixHandler& gh,
                                   unsigned long window_size_kb, unsigned long window_size);
        virtual void processVariant(const Variant& variant, const std::vector<Allele>& alleles);
        virtual bool isOfInterest(const Variant& position);
        void cleanup();
    };
}

#endif //SRC_VCF_PREDICTING_HANDLER_H
