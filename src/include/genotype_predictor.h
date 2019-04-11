#ifndef SRC_GENOTYPE_PREDICTOR_H
#define SRC_GENOTYPE_PREDICTOR_H

#include <vector>
#include "vcf_primitives.h"
#include <random>
#include <algorithm>
#include <exception>
#include <cmath>

namespace {
    typedef std::mt19937 Random;

    class Node {
        std::vector<double> class_weights;
    public:
        explicit Node(std::vector<double>&& class_weights) :class_weights(std::move(class_weights)){}

        double weight(size_t cl) {
            return class_weights[cl];
        }
    };

    class InnerNode : public Node {
        typedef std::unique_ptr<Node> nodePtr;
        nodePtr left_node;
        nodePtr right_node;
        vcf::AlleleType separator;
    public:
        InnerNode(std::vector<double>&& class_weights, nodePtr&& left_node, nodePtr&& right_node, vcf::AlleleType sep)
            :Node(std::move(class_weights)), left_node(std::move(left_node)), right_node(std::move(right_node)),
             separator(sep) {}
    };

    class LeafNode : public Node {
    public:
        explicit LeafNode(std::vector<double>&& class_weights) :Node(std::move(class_weights)) {};
    };

    int to_int(vcf::AlleleType type) {
        switch(type) {
            case vcf::HOMREF:
                return 0;
            case vcf::HET:
                return 1;
            case vcf::HOM:
                return 2;
            case vcf::MISSING:
                return 3;
        }
    }

    class Sample {
        int num;
        double w;
    public:
        Sample(int num, double weight) :num(num), w(weight){}

        int sample() {
            return num;
        }

        double weight() {
            return w;
        }
    };

    class Bags {
        std::vector<int> bags;
        std::vector<double> weights;
        double weights_sum;

    public:
        Bags(size_t size, Random& random) {
            for (int i = 0; i < size; i++) {
                bags.push_back(random() % size);
                weights.resize(size, 1.0);
            }
            weights_sum = size;
        }

        Bags() = default;

        void add(int sample, double weight) {
            bags.push_back(sample);
            weights.push_back(weight);
            weights_sum += weight;
        }

        std::vector<Sample> list() {
            std::vector<Sample> ret;
            ret.reserve(bags.size());
            for (int i = 0; i < bags.size(); i++) {
                ret.emplace_back(bags[i], weights[i]);
            }
            return ret;
        }

        double sum() {
            return weights_sum;
        }
    };
}

namespace vcf {
    typedef std::vector<std::vector<AlleleType>> Features;
    typedef std::vector<AlleleType> Labels;
    typedef std::pair<Bags, Bags> Split;

    class DecisionTree {
        const double EPS = 1e-8;
        Random random;
        std::unique_ptr<Node> root;
    public:
        DecisionTree(std::mt19937& random) :random(random) {}

        void fit(Features& features, Labels& y) {
        }

    private:
        std::tuple<double, double, double> counts(Bags& bags, Labels& labels) {
            double hom = 0.0;
            double het = 0.0;
            double alt = 0.0;
            for (auto s: bags.list()) {
                switch (labels[s.sample()]){
                    case HOM:
                        alt += s.weight();
                        break;
                    case HET:
                        het += s.weight();
                        break;
                    case HOMREF:
                        hom += s.weight();
                        break;
                }
            }
            return {hom, het, alt};
        }

        std::tuple<double, double, double> ratios(Bags& bags, Labels& labels) {
            auto count = counts(bags, labels);

            double hom = std::get<0>(count);
            double het = std::get<1>(count);
            double alt = std::get<2>(count);

            double sum = hom + het + alt;
            return {hom / sum, het / sum, alt / sum};
        }

        std::pair<Bags, Bags> split(Bags& curr, AlleleType splitBy, std::vector<AlleleType>& features, Labels& labels) {
            Bags left, right;
            auto list = curr.list();
            auto rs = ratios(curr, labels);
            double left_ratio = std::get<0>(rs);
            if (splitBy == HET) {
                left_ratio += std::get<1>(rs);
            }
            for (auto& el: list) {
                auto allele = features[el.sample()];
                if (allele == MISSING) {
                    left.add(el.sample(), el.weight() * left_ratio);
                    right.add(el.sample(), el.weight() * (1.0 - left_ratio));
                } else if (to_int(allele) <= to_int(splitBy)) {
                    left.add(el.sample(), el.weight());
                } else {
                    right.add(el.sample(), el.weight());
                }
            }
            return {left, right};
        }

        std::vector<int> sample(size_t n, size_t k, Random& random) {
            assert(k <= n);
            std::vector<int> arr(n);
            std::iota(arr.begin(), arr.end(), 0);
            std::shuffle(arr.begin(), arr.end(), random);
            return {arr.begin(), arr.begin() + k};
        }

        double score(Bags& samples, Labels& labels) {
            double info_gain = 0.0;
            auto cs_tuple = counts(samples, labels);
            std::vector<double> cs{std::get<0>(cs_tuple), std::get<1>(cs_tuple), std::get<2>(cs_tuple)};
            double sum = std::accumulate(cs.begin(), cs.end(), 0.0);
            for (double cnt: cs) {
                if (cnt != 0) {
                    double ratio = cnt / (double)sum;
                    info_gain -= ratio * std::log(ratio);
                }
            }
            return info_gain;
        }

        double split_score(Split& split, Labels& labels) {
            double sum = split.first.sum() + split.second.sum();
            double lr = split.first.sum() / sum;
            double rr = split.second.sum() / sum;
            return lr * score(split.first, labels) + rr * score(split.second, labels);
        }

        Node buildSubtree(Bags& bags, Features& features, Labels& values) {
            int k = std::floor(std::sqrt(features.size()));
            auto vars = sample(features.size(), k, random);
            int var_best = -1;
            AlleleType best_split = MISSING;
            double best_score = score(bags, values) - EPS;

            for (int var: vars) {
                auto hom_split = split(bags, HOMREF, features[var], values);
                auto het_split = split(bags, HET, features[var], values);
                double hom_score = split_score(hom_split, values);
                double het_score = split_score(het_split, values);
                if (hom_score < best_score) {
                    var_best = var;
                    best_split = HOMREF;
                    best_score = hom_score;
                }
                if (het_score < best_score) {
                    var_best = var;
                    best_split = HET;
                    best_score = het_score;
                }
            }

            if (best_split == MISSING) {
                auto cnts = counts(bags, values);
                std::vector<double> cs{std::get<0>(cnts), std::get<1>(cnts), std::get<2>(cnts)};
                return LeafNode(std::move(cs));
            }
        }
    };

}

#endif //SRC_GENOTYPE_PREDICTOR_H