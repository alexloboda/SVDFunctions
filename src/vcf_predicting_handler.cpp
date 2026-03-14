#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>

#include "include/vcf_predicting_handler.h"
#include "include/genotype_predictor.h"
#include "include/vcf_parser.h"

namespace {
    using std::size_t;

    static inline std::size_t safe_hardware_threads() {
        auto threads = (std::size_t)std::thread::hardware_concurrency();
        return std::max<std::size_t>(threads, 1);
    }

    static inline std::size_t metrics_threads_for(std::size_t tree_threads) {
        return std::max<std::size_t>(tree_threads / 4, 1);
    }

    static inline std::size_t pending_loo_limit_for(std::size_t sample_count, std::size_t metrics_threads) {
        if (sample_count >= 50000) {
            return 1;
        }
        return std::max<std::size_t>(metrics_threads, 1);
    }

    static inline int clamp_round_to_genotype(double x) {
        int r = (int)std::llround(x);
        if (r < 0) {
            return 0;
        }
        if (r > 2) {
            return 2;
        }
        return r;
    }

}

namespace vcf {
    PredictingHandler::PredictingHandler(const std::vector<std::string>& samples, GenotypeMatrixHandler& gh,
                                         int window_size_kb, int window_size,
                                         std::size_t rf_ntrees, unsigned int seed)
                                         :VariantsHandler(samples), curr_chr(-1), genotype_handler(gh), iterator{gh},
                                          window(window_size),
                                          thread_pool(safe_hardware_threads()),
                                          metrics_thread_pool(metrics_threads_for(safe_hardware_threads())),
                                          random_seed(seed),
                                          rf_ntrees(rf_ntrees),
                                          max_pending_loo_tasks(pending_loo_limit_for(samples.size(), metrics_threads_for(safe_hardware_threads()))) {
        auto variants = gh.desired_variants();
        int halfws = window_size_kb / 2;
        for (const Variant& v : variants) {
            Chromosome chr = v.position().chromosome();
            int pos = v.position().position();
            Range haplotype{chr, pos - halfws, pos + halfws};
            ranges.insert(haplotype);
        }
    }

    bool PredictingHandler::isOfInterest(const Variant& variant) {
        Position pos = variant.position();
        if (ranges.empty()) {
            return true;
        }
        return ranges.includes(pos);
    }

    void PredictingHandler::processVariant(const Variant& variant, std::shared_ptr<AlleleVector>& alleles,
                                           int64_t line_offset) {
        if (!isOfInterest(variant)) {
            return;
        }
        Position pos = variant.position();
        if (curr_chr.num() == -1) {
            curr_chr = pos.chromosome();
        }
        if (pos.chromosome() != curr_chr) {
            cleanup();
            curr_chr = pos.chromosome();
        }

        std::size_t logical_index = genotype_handler.next_row_position();
        if (genotype_handler.last_variant_consumed_row() && logical_index > 0) {
            --logical_index;
        }

        window.add(alleles, variant, line_offset, logical_index);
        if (window.is_full()) {
            while (iterator.dereferencable() && (*iterator).position().position() <= window.middle_point()) {
                auto dataset = window.dataset(*iterator);
                if (dataset.second.empty()) {
                    return;
                }
                Variant v = *iterator;
                fix_labels(v, std::move(dataset));
                ++iterator;
                collect_ready_loo_rows(false);
            }
        }
    }

    void PredictingHandler::cleanup() {
        for(; iterator.dereferencable(); ++iterator) {
            Variant var = *iterator;
            auto dataset = window.dataset(var);
            if (dataset.second.empty()) {
                break;
            }
            fix_labels(var, std::move(dataset));
            collect_ready_loo_rows(false);
        }
        collect_ready_loo_rows(true);
        window.clear();
    }

    std::size_t PredictingHandler::finalized_rows() const {
        return iterator.position();
    }

    int64_t PredictingHandler::resume_offset() const {
        return window.front_offset();
    }

    std::size_t PredictingHandler::resume_row_index() const {
        return window.front_logical_index();
    }

    void PredictingHandler::synchronize_loo() {
        collect_ready_loo_rows(true);
    }

    void PredictingHandler::discard_checkpointed_prefix(std::size_t new_row_offset) {
        iterator.discard_prefix(new_row_offset);
    }

    void PredictingHandler::discard_checkpointed_loo_prefix(std::size_t count) {
        if (count == 0) {
            return;
        }
        if (count > loo_rows.size()) {
            throw std::runtime_error("Attempted to discard more loo rows than are buffered");
        }
        loo_rows.erase(loo_rows.begin(), loo_rows.begin() + static_cast<std::ptrdiff_t>(count));
    }

    void PredictingHandler::collect_ready_loo_rows(bool wait_all) {
        while (!pending_loo_rows.empty()) {
            auto& next = pending_loo_rows.front();
            if (!wait_all && next.wait_for(std::chrono::seconds(0)) != std::future_status::ready) {
                break;
            }
            loo_rows.push_back(next.get());
            pending_loo_rows.pop_front();
        }
    }

    void PredictingHandler::collect_next_loo_row() {
        if (pending_loo_rows.empty()) {
            return;
        }
        loo_rows.push_back(pending_loo_rows.front().get());
        pending_loo_rows.pop_front();
    }

    void PredictingHandler::fix_labels(const Variant& variant, std::pair<Features, Labels> dataset) {
        const std::size_t logical_index = iterator.position();
        bool missing = false;
        for (auto l: dataset.second) {
            if (l == MISSING) {
                missing = true;
            }
        }
        if (!missing) {
            return;
        }
        TreeBuilder tree_builder = make_tree_builder(dataset);
        RandomForest forest{tree_builder, thread_pool, rf_ntrees, random_seed};

        // First pass: produce the returned genotype row (impute only missing).
        std::size_t n_observed = 0;
        std::size_t n_missing = 0;
        std::vector<size_t> observed_idx;
        observed_idx.reserve(dataset.second.size());

        std::vector<float> labels;
        labels.reserve(dataset.second.size());
        std::vector<AlleleType> sample_features(dataset.first.size());
        for (size_t i = 0; i < dataset.second.size(); i++) {
            AlleleType curr = dataset.second[i];
            if (curr == MISSING) {
                ++n_missing;
                for (size_t j = 0; j < dataset.first.size(); j++) {
                    sample_features[j] = dataset.first[j][i];
                }
                labels.push_back((float)forest.predict(sample_features));
            } else {
                ++n_observed;
                observed_idx.push_back(i);
                labels.push_back((float)to_int(curr));
            }
        }
        iterator.set(std::move(labels));

        std::size_t tree_pool_threads = 0;
        try {
            tree_pool_threads = thread_pool.n_threads();
        } catch (...) {
            tree_pool_threads = 0;
        }

        std::size_t metrics_pool_threads = 0;
        try {
            metrics_pool_threads = metrics_thread_pool.n_threads();
        } catch (...) {
            metrics_pool_threads = 0;
        }

        auto variant_name = (std::string)variant;
        pending_loo_rows.push_back(metrics_thread_pool.push([
            variant_name = std::move(variant_name),
            logical_index,
            n_observed,
            n_missing,
            tree_pool_threads,
            metrics_pool_threads,
            observed_idx = std::move(observed_idx),
            features = std::move(dataset.first),
            truth_labels = std::move(dataset.second),
            forest = std::move(forest)
        ]() mutable {
            struct EvalAccum {
                double sum_abs = 0.0;
                double sum_sq = 0.0;
                std::size_t correct = 0;
                std::size_t count = 0;
            };

            ImputationLooRow row;
            row.variant = std::move(variant_name);
            row.logical_index = logical_index;
            row.n_observed = n_observed;
            row.n_missing = n_missing;

            auto eval_rf_range = [&](size_t begin, size_t end) -> EvalAccum {
                EvalAccum acc;
                std::vector<AlleleType> sample_features(features.size());
                for (size_t k = begin; k < end; k++) {
                    size_t sample_index = observed_idx[k];
                    for (size_t j = 0; j < features.size(); j++) {
                        sample_features[j] = features[j][sample_index];
                    }
                    double pred = forest.predict_oob(sample_features, sample_index);
                    double truth = (double)to_int(truth_labels[sample_index]);
                    double err = pred - truth;
                    acc.sum_abs += std::abs(err);
                    acc.sum_sq += err * err;
                    int rounded = clamp_round_to_genotype(pred);
                    if ((double)rounded == truth) {
                        acc.correct++;
                    }
                    acc.count++;
                }
                return acc;
            };

            EvalAccum rf_acc;
            if (observed_idx.size() >= 512 && metrics_pool_threads >= 2) {
                const size_t target_tasks = std::min<std::size_t>(metrics_pool_threads, 8);
                const size_t chunk = std::max<std::size_t>((observed_idx.size() + target_tasks - 1) / target_tasks, 128);
                std::vector<std::future<EvalAccum>> futures;
                futures.reserve((observed_idx.size() + chunk - 1) / chunk);
                for (size_t begin = 0; begin < observed_idx.size(); begin += chunk) {
                    size_t end = std::min(begin + chunk, observed_idx.size());
                    futures.push_back(std::async(std::launch::async, [begin, end, &eval_rf_range]() {
                        return eval_rf_range(begin, end);
                    }));
                }
                for (auto& future: futures) {
                    auto acc = future.get();
                    rf_acc.sum_abs += acc.sum_abs;
                    rf_acc.sum_sq += acc.sum_sq;
                    rf_acc.correct += acc.correct;
                    rf_acc.count += acc.count;
                }
            } else {
                rf_acc = eval_rf_range(0, observed_idx.size());
            }

            if (rf_acc.count == 0) {
                row.oob_mae = std::numeric_limits<double>::quiet_NaN();
                row.oob_rmse = std::numeric_limits<double>::quiet_NaN();
                row.rounded_acc = std::numeric_limits<double>::quiet_NaN();
            } else {
                row.oob_mae = rf_acc.sum_abs / (double)rf_acc.count;
                row.oob_rmse = std::sqrt(rf_acc.sum_sq / (double)rf_acc.count);
                row.rounded_acc = (double)rf_acc.correct / (double)rf_acc.count;
            }

            return row;
        }));

        while (pending_loo_rows.size() > max_pending_loo_tasks) {
            collect_next_loo_row();
        }
    }

    TreeBuilder PredictingHandler::make_tree_builder(const std::pair<Features, Labels>& dataset) {
        size_t mtry = ceil(sqrt(dataset.first.size()));
        return {dataset.first, dataset.second, mtry};
    }

    Window::Window(size_t max_size) :max_size(max_size), start(0) {}

    void Window::clear() {
        features.clear();
        variants.clear();
        line_offsets.clear();
        logical_indices.clear();
        start = 0;
    }

    std::pair<Features, Labels> Window::dataset(const Variant& v) {
        Features fs;
        Labels lbls;
        size_t none = std::numeric_limits<size_t>::max();
        size_t curr_num = none;
        bool flipped = false;

        for (size_t i = 0; i < variants.size(); i++) {
            if (variants[i] == v || variants[i] == v.reversed()) {
                curr_num = i;
                if (variants[i] == v.reversed()) {
                    flipped = true;
                }
            }
        }
        if (curr_num == none) {
            return {};
        }
        lbls = features[curr_num]->vector(flipped);
        for (size_t i = 0; i < features.size(); i++) {
            if (i != curr_num) {
                try {
                    std::vector<AlleleType> row = features[i]->vector(flipped);
                    fs.push_back(std::move(row));
                } catch (ParserException& e) {}
            }
        }
        return {std::move(fs), std::move(lbls)};
    }

    void Window::add(std::shared_ptr<AlleleVector>& alleles, const Variant& variant,
                     int64_t line_offset, std::size_t logical_index) {
        if (features.size() < max_size) {
            variants.push_back(variant);
            features.push_back(alleles);
            line_offsets.push_back(line_offset);
            logical_indices.push_back(logical_index);
        } else {
            variants.pop_front();
            features.pop_front();
            line_offsets.pop_front();
            logical_indices.pop_front();
            variants.push_back(variant);
            features.push_back(alleles);
            line_offsets.push_back(line_offset);
            logical_indices.push_back(logical_index);
        }
    }

    int Window::middle_point() {
        if (variants.empty()) {
            return -1;
        }
        return variants[variants.size() / 2].position().position();
    }

    bool Window::is_full() {
        return features.size() == max_size;
    }

    bool Window::empty() const {
        return variants.empty();
    }

    int64_t Window::front_offset() const {
        if (line_offsets.empty()) {
            return -1;
        }
        return line_offsets.front();
    }

    std::size_t Window::front_logical_index() const {
        if (logical_indices.empty()) {
            return 0;
        }
        return logical_indices.front();
    }
}
