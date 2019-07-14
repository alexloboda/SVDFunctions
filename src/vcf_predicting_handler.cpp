#include "include/vcf_predicting_handler.h"
#include "include/genotype_predictor.h"
#include "include/vcf_parser.h"

namespace {
    using std::size_t;
}

namespace vcf {
    PredictingHandler::PredictingHandler(const std::vector<std::string>& samples, GenotypeMatrixHandler& gh,
                                         int window_size_kb, int window_size)
                                         :VariantsHandler(samples), curr_chr(-1), iterator{gh},
                                          window(window_size),
                                          thread_pool(std::thread::hardware_concurrency()) {
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
                fix_labels(dataset);
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
            fix_labels(dataset);
        }
        window.clear();
    }

    void PredictingHandler::fix_labels(const std::pair<Features, Labels>& dataset) {
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
        for (size_t i = 0; i < variants.size(); i++) {
            if (variants[i] == v) {
                curr_num = i;
            }
        }
        if (curr_num == none) {
            return {};
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
