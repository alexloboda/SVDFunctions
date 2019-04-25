#ifndef SRC_VCF_PREDICTING_HANDLER_H
#define SRC_VCF_PREDICTING_HANDLER_H

#include "vcf_handlers.h"
#include "genotype_predictor.h"
#include "vcf_parser.h"

namespace vcf {
    class Window {
        std::deque<std::shared_ptr<AlleleVector>> features;
        std::deque<Variant> variants;
        size_t max_size;
        size_t start;
    public:
        explicit Window(size_t max_size);
        void clear();
        void add(std::shared_ptr<AlleleVector>& alleles, const Variant& variant);
        std::pair<Features, Labels> dataset(const Variant& v);
        int middle_point();
        bool is_full();
    };

    class PredictingHandler : public VariantsHandler {
        Chromosome curr_chr;
        std::unordered_map<int, std::set<Range>> ranges;
        GenotypeMatrixIterator iterator;
        Window window;
        std::vector<size_t> dataset_samples;
        cxxpool::thread_pool thread_pool;

        TreeBuilder make_tree_builder(const std::pair<Features, Labels>& dataset);
    public:
        explicit PredictingHandler(const std::vector<std::string>& samples, GenotypeMatrixHandler& gh,
                                   int window_size_kb, int window_size, std::vector<size_t>& dataset_samples);
        void processVariant(const Variant& variant, std::shared_ptr<AlleleVector>& alleles) override;
        bool isOfInterest(const Variant& position) override;
        void cleanup();
        void fix_labels(std::pair<Features, Labels>& dataset);
    };
}

#endif //SRC_VCF_PREDICTING_HANDLER_H
