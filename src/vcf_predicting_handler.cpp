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
                                         int window_size_kb, int window_size)
                                         :VariantsHandler(samples), gh(gh), curr_chr(-1), iterator{gh},
                                          window(window_size, window_size_kb){
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

        std::vector<AlleleType> sample;
        std::for_each(alleles->begin(), alleles->end(), [&sample](const Allele& allele){
            sample.push_back(allele.alleleType());
        });
        window.add(sample, variant);
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
        size_t mtry = ceil(sqrt(dataset.first.size()));
        TreeBuilder tree_builder{dataset.first, dataset.second, mtry};
        RandomForest forest{tree_builder};
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

    Window::Window(size_t max_size, size_t max_size_kb) :max_size(max_size), max_size_kb(max_size_kb), start(0) {}

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
        lbls = features[curr_num];
        for (size_t i = 0; i < features.size(); i++) {
            if (i != curr_num) {
                std::vector<AlleleType> row;
                fs.push_back(features[i]);
            }
        }
        return {std::move(fs), std::move(lbls)};
    }

    void Window::add(const std::vector<AlleleType>& alleles, const Variant& variant) {
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
