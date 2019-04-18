#ifndef SRC_VCF_PREDICTING_HANDLER_H
#define SRC_VCF_PREDICTING_HANDLER_H

#include "vcf_handlers.h"
#include "genotype_predictor.h"

namespace vcf {
    class Window {
        std::deque<std::vector<AlleleType>> features;
        std::deque<Variant> variants;
        size_t max_size;
        size_t max_size_kb;
        size_t start;
    public:
        explicit Window(size_t max_size, size_t max_size_kb);
        void clear();
        void add(const std::vector<AlleleType>& alleles, const Variant& variant);
        std::pair<Features, Labels> dataset(const Variant& v);
        int middle_point();
        bool is_full();
    };

    class PredictingHandler : public VariantsHandler {
        GenotypeMatrixHandler& gh;
        Chromosome curr_chr;
        std::unordered_map<int, std::set<Range>> ranges;
        GenotypeMatrixIterator iterator;
        Window window;

    public:
        explicit PredictingHandler(const std::vector<std::string>& samples, GenotypeMatrixHandler& gh,
                                   int window_size_kb, int window_size);
        void processVariant(const Variant& variant, const std::vector<Allele>& alleles) override;
        bool isOfInterest(const Variant& position) override;
        void cleanup();
        void fix_labels(std::pair<Features, Labels>& dataset);
    };
}

#endif //SRC_VCF_PREDICTING_HANDLER_H
