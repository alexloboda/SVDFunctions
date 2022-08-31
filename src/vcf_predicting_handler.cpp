#include "include/vcf_predicting_handler.h"
#include "include/genotype_predictor.h"
#include "include/vcf_parser.h"

#include <random>
#include <utility>

namespace {
    using std::size_t;
}

namespace vcf {
    PredictingHandler::PredictingHandler(const std::vector<std::string>& samples, GenotypeMatrixHandler& gh,
                                         int window_size_kb, int window_size, bool cv, size_t trees, bool mean)
                                         :VariantsHandler(samples), curr_chr(-1), iterator{gh},
                                          window(window_size),
                                          thread_pool(std::thread::hardware_concurrency()), cv(cv),
                                          trees(trees), mean(mean) {
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

    namespace {
        std::vector<AlleleType> subset(const std::vector<AlleleType>& set, const std::vector<size_t>& idx) {
            std::vector<AlleleType> ret;
            for (auto i: idx) {
                ret.push_back(set[i]);
            }
            return ret;
        }

        std::pair<Features, Labels> subdataset(std::pair<Features, Labels> dataset, const std::vector<size_t>& idx) {
            dataset.second = subset(dataset.second, idx);
            for (size_t i = 0; i < dataset.first.size(); i++) {
                dataset.first[i] = subset(dataset.first[i], idx);
            }
            return dataset;
        }
    }

    void PredictingHandler::cross_validation(const std::pair<Features, Labels>& dataset) {
        std::vector<size_t> nonmissing;
        for (size_t i = 0; i < dataset.second.size(); i++) {
            if (dataset.second[i] != MISSING) {
                nonmissing.push_back(i);
            }
        }
        int seed = rand();
        Random random(seed);
        std::shuffle(nonmissing.begin(), nonmissing.end(), random);
        int chunk = nonmissing.size() / 5;
        for (int k = 0; k < 5; k++) {
            std::vector<double> mse = {0.0, 0.0, 0.0, 0.0};
            std::vector<size_t> train;
            std::vector<size_t> test;
            for (size_t i = 0; i < nonmissing.size(); i++) {
                if (i >= chunk * k && i < chunk * (k + 1)) {
                    test.push_back(nonmissing[i]);
                } else {
                    train.push_back(nonmissing[i]);
                }
            }
            auto training_set = subdataset(dataset, train);

            TreeBuilder tree_builder = make_tree_builder(training_set, mean);
            RandomForest forest{tree_builder, thread_pool, trees, !mean};
            std::vector<int> counts = {0, 0, 0};
            for (size_t i: test) {
                std::vector<AlleleType> features;
                for (size_t j = 0; j < dataset.first.size(); j++) {
                    features.push_back(dataset.first[j][i]);
                }
                auto predicted = forest.predict(features);
                double error = predicted - to_int(dataset.second[i]);
                mse[to_int(dataset.second[i])] += error * error;
                counts[to_int(dataset.second[i])]++;
            }
            mse[3] = (mse[0] + mse[1] + mse[2]) / (double)test.size();
            for (int i = 0; i < 3; i++) {
                mse[i] /= counts[i];
            }
            cv_mse.push_back(mse);
        }
    }

    void PredictingHandler::fix_labels(const std::pair<Features, Labels>& dataset) {
        if (cv) {
            cross_validation(dataset);
        }
        bool missing = false;
        for (auto l: dataset.second) {
            if (l == MISSING) {
                missing = true;
            }
        }
        if (!missing) {
            return;
        }
        TreeBuilder tree_builder = make_tree_builder(dataset, false);
        RandomForest forest{tree_builder, thread_pool, trees};
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

    TreeBuilder PredictingHandler::make_tree_builder(const std::pair<Features, Labels>& dataset, bool mean) {
        size_t mtry = ceil(sqrt(dataset.first.size()));
        if (mean) {
            return {{}, dataset.second, mtry};
        } else {
            return {dataset.first, dataset.second, mtry};
        }
    }

std::vector<std::vector<double>> PredictingHandler::mses() const {
    return cv_mse;
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
                try {
                    std::vector<AlleleType> row = features[i]->vector();
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
