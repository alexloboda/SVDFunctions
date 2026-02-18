#include <cmath>
#include <limits>
#include <algorithm>
#include <future>

#include "include/vcf_predicting_handler.h"
#include "include/genotype_predictor.h"
#include "include/vcf_parser.h"

namespace {
    using std::size_t;

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
                                         int window_size_kb, int window_size, unsigned int seed)
                                         :VariantsHandler(samples), curr_chr(-1), iterator{gh},
                                          window(window_size),
                                          thread_pool(std::thread::hardware_concurrency()),
                                          random_seed(seed) {
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

    void PredictingHandler::processVariant(const Variant& variant, std::shared_ptr<AlleleVector>& alleles) {
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

        window.add(alleles, variant);
        if (window.is_full()) {
            while (iterator.dereferencable() && (*iterator).position().position() <= window.middle_point()) {
                auto dataset = window.dataset(*iterator);
                if (dataset.second.empty()) {
                    return;
                }
                Variant v = *iterator;
                fix_labels(v, dataset);
                ++iterator;
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
            fix_labels(var, dataset);
        }
        window.clear();
    }

    void PredictingHandler::fix_labels(const Variant& variant, const std::pair<Features, Labels>& dataset) {
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
        RandomForest forest{tree_builder, thread_pool, /*ntrees=*/100, random_seed};

        // First pass: produce the returned genotype row (impute only missing).
        std::size_t n_observed = 0;
        std::size_t n_missing = 0;
        std::vector<size_t> observed_idx;
        observed_idx.reserve(dataset.second.size());

        std::vector<float> labels;
        labels.reserve(dataset.second.size());
        for (size_t i = 0; i < dataset.second.size(); i++) {
            AlleleType curr = dataset.second[i];
            if (curr == MISSING) {
                ++n_missing;
                std::vector<AlleleType> features;
                features.reserve(dataset.first.size());
                for (size_t j = 0; j < dataset.first.size(); j++) {
                    features.push_back(dataset.first[j][i]);
                }
                labels.push_back((float)forest.predict(features));
            } else {
                ++n_observed;
                observed_idx.push_back(i);
                labels.push_back((float)to_int(curr));
            }
        }
        iterator.set(labels);

        // RF OOB evaluation on observed genotypes (parallelized in chunks when beneficial).
        struct EvalAccum {
            double sum_abs = 0.0;
            double sum_sq = 0.0;
            std::size_t correct = 0;
            std::size_t count = 0;
        };

        auto eval_rf_range = [&](size_t begin, size_t end) -> EvalAccum {
            EvalAccum acc;
            for (size_t k = begin; k < end; k++) {
                size_t i = observed_idx[k];
                std::vector<AlleleType> features;
                features.reserve(dataset.first.size());
                for (size_t j = 0; j < dataset.first.size(); j++) {
                    features.push_back(dataset.first[j][i]);
                }
                double pred = forest.predict_oob(features, i);
                double truth = (double)to_int(dataset.second[i]);
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

        std::size_t pool_threads = 0;
        try {
            pool_threads = thread_pool.n_threads();
        } catch (...) {
            pool_threads = 0;
        }

        EvalAccum rf_acc;
        if (observed_idx.size() >= 512 && pool_threads >= 2) {
            const size_t target_tasks = std::min<std::size_t>(pool_threads * 4, 32);
            const size_t chunk = std::max<std::size_t>((observed_idx.size() + target_tasks - 1) / target_tasks, 128);
            std::vector<std::future<EvalAccum>> futures;
            for (size_t begin = 0; begin < observed_idx.size(); begin += chunk) {
                size_t end = std::min(begin + chunk, observed_idx.size());
                futures.push_back(thread_pool.push([&eval_rf_range, begin, end]() { return eval_rf_range(begin, end); }));
            }
            for (auto& f : futures) {
                auto a = f.get();
                rf_acc.sum_abs += a.sum_abs;
                rf_acc.sum_sq += a.sum_sq;
                rf_acc.correct += a.correct;
                rf_acc.count += a.count;
            }
        } else {
            rf_acc = eval_rf_range(0, observed_idx.size());
        }


        ImputationLooRow row;
        row.variant = (std::string)variant;
        row.n_observed = n_observed;
        row.n_missing = n_missing;
        if (rf_acc.count == 0) {
            row.oob_mae = std::numeric_limits<double>::quiet_NaN();
            row.oob_rmse = std::numeric_limits<double>::quiet_NaN();
            row.rounded_acc = std::numeric_limits<double>::quiet_NaN();
        } else {
            row.oob_mae = rf_acc.sum_abs / (double)rf_acc.count;
            row.oob_rmse = std::sqrt(rf_acc.sum_sq / (double)rf_acc.count);
            row.rounded_acc = (double)rf_acc.correct / (double)rf_acc.count;
        }

        loo_rows.push_back(std::move(row));
    }

    TreeBuilder PredictingHandler::make_tree_builder(const std::pair<Features, Labels>& dataset) {
        size_t mtry = ceil(sqrt(dataset.first.size()));
        return {dataset.first, dataset.second, mtry};
    }

    Window::Window(size_t max_size) :max_size(max_size), start(0) {}

    void Window::clear() {
        features.clear();
        variants.clear();
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

    void Window::add(std::shared_ptr<AlleleVector>& alleles, const Variant& variant) {
        if (features.size() < max_size) {
            variants.push_back(variant);
            features.push_back(alleles);
        } else {
            variants.pop_front();
            features.pop_front();
            variants.push_back(variant);
            features.push_back(alleles);
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
}
