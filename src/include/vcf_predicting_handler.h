#ifndef SRC_VCF_PREDICTING_HANDLER_H
#define SRC_VCF_PREDICTING_HANDLER_H

#include <deque>
#include <future>
#include <string>

#include "vcf_handlers.h"
#include "genotype_predictor.h"
#include "vcf_parser.h"

namespace vcf {
    struct ImputationLooRow {
        std::string variant;
        std::size_t logical_index;
        std::size_t n_observed;
        std::size_t n_missing;
        double oob_mae;
        double oob_rmse;
        double rounded_acc;
    };

    class Window {
        std::deque<std::shared_ptr<AlleleVector>> features;
        std::deque<Variant> variants;
        std::deque<int64_t> line_offsets;
        std::deque<std::size_t> logical_indices;
        std::size_t max_size;
        std::size_t start;
    public:
        explicit Window(std::size_t max_size);
        void clear();
        void add(std::shared_ptr<AlleleVector>& alleles, const Variant& variant,
                 int64_t line_offset, std::size_t logical_index);
        std::pair<Features, Labels> dataset(const Variant& v);
        int middle_point();
        bool is_full();
        bool empty() const;
        int64_t front_offset() const;
        std::size_t front_logical_index() const;
    };

    class PredictingHandler : public VariantsHandler {
        Chromosome curr_chr;
        RangeSet ranges;
        GenotypeMatrixHandler& genotype_handler;
        GenotypeMatrixIterator iterator;
        Window window;
        cxxpool::thread_pool thread_pool;
        cxxpool::thread_pool metrics_thread_pool;
        unsigned int random_seed;
        std::size_t rf_ntrees;
        std::size_t max_pending_loo_tasks;

        std::deque<std::future<ImputationLooRow>> pending_loo_rows;
        std::vector<ImputationLooRow> loo_rows;

        TreeBuilder make_tree_builder(const std::pair<Features, Labels>& dataset);
        void collect_ready_loo_rows(bool wait_all);
        void collect_next_loo_row();
    public:
        explicit PredictingHandler(const std::vector<std::string>& samples, GenotypeMatrixHandler& gh,
                                   int window_size_kb, int window_size,
                                   std::size_t rf_ntrees = 50,
                                   unsigned int seed = 42);
        void processVariant(const Variant& variant, std::shared_ptr<AlleleVector>& alleles,
                            int64_t line_offset = -1) override;
        bool isOfInterest(const Variant& position) override;
        void cleanup();

        void fix_labels(const Variant& variant, std::pair<Features, Labels> dataset);

        std::size_t finalized_rows() const;
        int64_t resume_offset() const;
        std::size_t resume_row_index() const;
        void synchronize_loo();
        void discard_checkpointed_prefix(std::size_t new_row_offset);
        void discard_checkpointed_loo_prefix(std::size_t count);

        const std::vector<ImputationLooRow>& imputation_loo() const { return loo_rows; }
    };
}

#endif //SRC_VCF_PREDICTING_HANDLER_H
