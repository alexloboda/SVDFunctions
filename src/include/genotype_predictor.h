#ifndef SRC_GENOTYPE_PREDICTOR_H
#define SRC_GENOTYPE_PREDICTOR_H

#include <vector>
#include <random>
#include <algorithm>
#include <exception>
#include <cmath>

#include "vcf_primitives.h"
#include "third-party/cxxpool.h"

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
        virtual double predict(std::vector<vcf::AlleleType>& features) const = 0;
        virtual ~Node() = default;

        std::vector<double> weights();

    protected:
        static double prediction(const std::vector<double>& alpha);
    };

    typedef std::shared_ptr<Node> NodePtr;

    class DecisionTree {
        NodePtr root;
    public:
        DecisionTree(const DecisionTree&) = delete;
        DecisionTree(DecisionTree&& other) noexcept;
        static constexpr double EPS = 1e-8;

        double predict(std::vector<AlleleType>& features) const;
        double accuracy();

        friend class TreeBuilder;
    private:
        explicit DecisionTree(NodePtr root);
    };

    class Bags;

    class TreeBuilder {
        const Features features;
        const Labels values;
        std::size_t max_features;
    public:
        TreeBuilder(const Features&, const Labels& labels, std::size_t max_features);
        DecisionTree build_a_tree(Random& random, bool bagging = true) const;
    private:
        NodePtr buildSubtree(const Bags& bags, Random& random) const;

    };

    class RandomForest {
        std::vector<DecisionTree> predictors;
    public:
        RandomForest(const TreeBuilder& treeBuilder, cxxpool::thread_pool& thread_pool, std::size_t trees = 12, unsigned int seed = 42);
        double predict(std::vector<AlleleType>& features);
    };
}

#endif //SRC_GENOTYPE_PREDICTOR_H