#ifndef SRC_VCF_PREDICTING_HANDLER_H
#define SRC_VCF_PREDICTING_HANDLER_H

#include "vcf_handlers.h"
#include "genotype_predictor.h"

namespace vcf {
    class Window {
        Features features;
        std::vector<Variant> variants;
        size_t max_size;
    public:
        explicit Window(size_t max_size);
        void clear();
        void add(std::vector<AlleleType>& alleles, Variant& variant);
        std::pair<Features, Labels> dataset(const Variant& v);
    };

    class PredictingHandler : public VariantsHandler {
        const int window_size_kb;
        const int window_size;

        GenotypeMatrixHandler& gh;
        Chromosome curr_chr;
        std::unordered_map<int, std::set<Range>> ranges;
        GenotypeMatrixIterator iterator;

    public:
        explicit PredictingHandler(const std::vector<std::string>& samples, GenotypeMatrixHandler& gh,
                                   int window_size_kb, int window_size);
        void processVariant(const Variant& variant, const std::vector<Allele>& alleles) override;
        bool isOfInterest(const Variant& position) override;
        void cleanup();
    };
}

#endif //SRC_VCF_PREDICTING_HANDLER_H
