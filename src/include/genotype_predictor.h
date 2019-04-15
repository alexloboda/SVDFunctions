#ifndef SRC_GENOTYPE_PREDICTOR_H
#define SRC_GENOTYPE_PREDICTOR_H

#include <vector>
#include "vcf_primitives.h"
#include <random>
#include <algorithm>
#include <exception>
#include <cmath>

namespace vcf {
    typedef std::vector<std::vector<vcf::AlleleType>> Features;
    typedef std::vector<vcf::AlleleType> Labels;
    typedef std::mt19937 Random;

    class Node {
    protected:
        std::vector<double> class_weights;
        double acc;

    public:
        explicit Node(std::vector<double>&& class_weights);

        double accuracy();
        virtual double predict(std::vector<vcf::AlleleType>& features) = 0;
        virtual ~Node() = default;

        std::vector<double> weights();

    protected:
        static double prediction(std::vector<double>& alpha);
    };

    typedef std::shared_ptr<Node> NodePtr;

    class DecisionTree {
        NodePtr root;
    public:
        DecisionTree(const DecisionTree&) = delete;
        DecisionTree(DecisionTree&& other) noexcept;
        static constexpr double EPS = 1e-8;

        double predict(std::vector<AlleleType>& features);
        double accuracy();

        friend class TreeBuilder;
    private:
        explicit DecisionTree(NodePtr root);
    };

    class Bags;

    class TreeBuilder {
        const Features features;
        const Labels values;
        size_t max_features;
    public:
        TreeBuilder(Features& features, Labels& labels, size_t max_features);
        DecisionTree build_a_tree(Random& random);
    private:
        NodePtr buildSubtree(Bags& bags, Random& random);

    };

}

#endif //SRC_GENOTYPE_PREDICTOR_H