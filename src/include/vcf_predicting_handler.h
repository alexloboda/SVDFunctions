#ifndef SRC_VCF_PREDICTING_HANDLER_H
#define SRC_VCF_PREDICTING_HANDLER_H

#include <string>

#include "vcf_handlers.h"
#include "genotype_predictor.h"
#include "vcf_parser.h"

namespace vcf {
    struct ImputationLooRow {
        std::string variant;
        std::size_t n_observed;
        std::size_t n_missing;
        double oob_mae;
        double oob_rmse;
        double rounded_acc;
    };

    class Window {
        std::deque<std::shared_ptr<AlleleVector>> features;
        std::deque<Variant> variants;
        std::size_t max_size;
        std::size_t start;
    public:
        explicit Window(std::size_t max_size);
        void clear();
        void add(std::shared_ptr<AlleleVector>& alleles, const Variant& variant);
        std::pair<Features, Labels> dataset(const Variant& v);
        int middle_point();
        bool is_full();
    };

    class PredictingHandler : public VariantsHandler {
        Chromosome curr_chr;
        RangeSet ranges;
        GenotypeMatrixIterator iterator;
        Window window;
        cxxpool::thread_pool thread_pool;
        unsigned int random_seed; 

        std::vector<ImputationLooRow> loo_rows;

        TreeBuilder make_tree_builder(const std::pair<Features, Labels>& dataset);
    public:
        explicit PredictingHandler(const std::vector<std::string>& samples, GenotypeMatrixHandler& gh,
                                   int window_size_kb, int window_size, unsigned int seed = 42); 
        void processVariant(const Variant& variant, std::shared_ptr<AlleleVector>& alleles) override;
        bool isOfInterest(const Variant& position) override;
        void cleanup();

        void fix_labels(const Variant& variant, const std::pair<Features, Labels>& dataset);

        const std::vector<ImputationLooRow>& imputation_loo() const { return loo_rows; }
    };
}

#endif //SRC_VCF_PREDICTING_HANDLER_H
