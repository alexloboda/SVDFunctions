#include "include/vcf_predicting_handler.h"
#include "include/genotype_predictor.h"
#include "include/vcf_parser.h"

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
                                         int window_size_kb, int window_size, std::vector<size_t>& dataset_samples)
                                         :VariantsHandler(samples), curr_chr(-1), iterator{gh},
                                          window(window_size, window_size_kb), dataset_samples(dataset_samples),
                                          thread_pool(std::thread::hardware_concurrency()) {
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
        auto& set = ranges[pos.chromosome().num()];
        // Workaround, c++11
        auto it = set.lower_bound(Range(pos.chromosome(), -1, pos.position()));
        return it != set.end() && it->includes(pos);
    }

    void PredictingHandler::processVariant(const Variant& variant, std::shared_ptr<AlleleVector>& alleles) {
        if (!isOfInterest(variant)) {
            return;
        }
        Position pos = variant.position();
        if (pos.chromosome() != curr_chr) {
            cleanup();
            curr_chr = pos.chromosome();
        }

        window.add(alleles, variant);
        if (window.is_full()) {
            while (iterator.dereferencable() && (*iterator).position().position() <= window.middle_point()) {
                auto dataset = window.dataset(*iterator);
                fix_labels(dataset);
                ++iterator;
            }
        }

    }

    void PredictingHandler::cleanup() {
        for(; iterator.dereferencable(); ++iterator) {
            Variant var = *iterator;
            auto dataset = window.dataset(var);
            fix_labels(dataset);
        }
        window.clear();
    }

    void PredictingHandler::fix_labels(std::pair<Features, Labels>& dataset) {
        TreeBuilder tree_builder = make_tree_builder(dataset);
        RandomForest forest{tree_builder, thread_pool};
        std::vector<float> labels;
        for (size_t i = 0; i < dataset.second.size(); i++) {
            AlleleType curr = dataset.second[i];
            if (curr == MISSING) {
                std::vector<AlleleType> features;
                for (size_t j = 0; j < dataset.first.size(); j++) {
                    features.push_back(dataset.first[j][i]);
                }
                labels.push_back(forest.predict(features));
            } else {
                labels.push_back(to_int(curr));
            }
        }
        iterator.set(labels);
    }

    TreeBuilder PredictingHandler::make_tree_builder(std::pair<Features, Labels>& dataset) {
        if (dataset_samples.size() == dataset.second.size()) {
            size_t mtry = ceil(sqrt(dataset.first.size()));
            return {std::move(dataset.first), std::move(dataset.second), mtry};
        } else {
            Features features;
            Labels labels;
            for (size_t i = 0; i < features.size(); i++) {
                features.emplace_back();
                for (size_t s: dataset_samples) {
                    features[i].push_back(dataset.first[i][s]);
                }
            }
            for (size_t s: dataset_samples) {
                labels.push_back(labels[s]);
            }
            size_t mtry = ceil(sqrt(features.size()));
            return {std::move(features), std::move(labels), mtry};
        }
    }

    Window::Window(size_t max_size, size_t max_size_kb) :max_size(max_size), start(0) {}

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
        for (size_t i = 0; i < variants.size(); i++) {
            if (variants[i] == v) {
                curr_num = i;
            }
        }
        if (curr_num == none) {
            throw std::logic_error("No values for training set. Potentially unreachable code.");
        }
        lbls = features[curr_num]->vector();
        for (size_t i = 0; i < features.size(); i++) {
            if (i != curr_num) {
                std::vector<AlleleType> row;
                fs.push_back(features[i]->vector());
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
