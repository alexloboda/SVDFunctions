#ifndef SRC_GENOTYPE_PREDICTOR_H
#define SRC_GENOTYPE_PREDICTOR_H

#include <vector>
#include "vcf_primitives.h"
#include <random>
#include <algorithm>
#include <exception>

namespace {
    typedef std::mt19937 Random;
    const size_t CLASSES = 3;

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

    class LeafNode : public Node {};

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

    public:
        Bags(size_t size, Random& random) {
            for (int i = 0; i < size; i++) {
                bags.push_back(random() % size);
                weights.resize(size, 1.0);
            }
        }

        Bags() = default;

        void add(int sample, double weight) {
            bags.push_back(sample);
            weights.push_back(weight);
        }

        std::vector<Sample> list() {
            std::vector<Sample> ret;
            ret.reserve(bags.size());
            for (int i = 0; i < bags.size(); i++) {
                ret.emplace_back(bags[i], weights[i]);
            }
            return ret;
        }
    };
}

namespace vcf {
    typedef std::vector<std::vector<AlleleType>> Features;
    typedef std::vector<AlleleType> Labels;

    class DecisionTree {
        Random random;
        std::unique_ptr<Node> root;
    public:
        DecisionTree(std::mt19937& random) :random(random) {}

        void fit(Features& features, Labels& y) {
        }

    private:
        Bags split(int var, Bags& curr, AlleleType splitBy, Features& features, Labels& labels) {
            Bags left, right;
            auto list = curr.list();
            auto& var_data = features[var];
            for (auto& el: list) {

            }
        }

        std::vector<int> sample(size_t n, size_t k, Random& random) {
            assert(k <= n);
            std::vector<int> arr(n);
            std::iota(arr.begin(), arr.end(), 0);
            std::shuffle(arr.begin(), arr.end(), random);
            return {arr.begin(), arr.begin() + k};
        }

        Node buildSubtree(Bags& bags, Features& features, Labels& values) {
            int k = std::floor(std::sqrt(features.size()));
            auto vars = sample(features.size(), k, random);
            for (int var: vars) {

            }
        }
    };

}

#endif //SRC_GENOTYPE_PREDICTOR_H