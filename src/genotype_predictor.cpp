#include <utility>
#include <random>
#include <iostream>

#include "include/genotype_predictor.h"
#include "include/vcf_primitives.h"
#include "include/third-party/cxxpool.h"

namespace {
    using std::size_t;

    class Sample {
        int num;
        double w;
    public:
        Sample(int num, double weight) :num(num), w(weight){}

        int sample() const {
            return num;
        }

        double weight() const {
            return w;
        }
    };

    const double EPS = 1e-8;
}

namespace vcf {

    class Bags {
        std::vector<Sample> samples;
        double weights_sum = 0.0;

    public:
        Bags(const Bags& bags, Random& random) {
            std::vector<double> prefix_weights;
            double curr = 0.0;
            auto list = bags.list();
            for (auto s: list) {
                curr += s.weight();
                prefix_weights.push_back(curr);
            }
            // just to be sure lower_bound won't return end()
            prefix_weights[prefix_weights.size() - 1] += 1.0;
            std::uniform_real_distribution<double> r(0.0, curr);
            for (size_t i = 0; i < list.size(); i++) {
                auto it = std::lower_bound(prefix_weights.begin(), prefix_weights.end(), r(random));
                size_t pos = std::distance(prefix_weights.begin(), it);

                Sample insertion = list[pos];
                samples.push_back(insertion);
                weights_sum += insertion.weight();
            }
        }

        Bags() = default;
        Bags(const Bags&) = delete;
        Bags& operator=(const Bags&) = delete;
        Bags(Bags&& other) noexcept :samples(std::move(other.samples)), weights_sum(other.weights_sum) {}

        void add(int sample, double weight) {
            samples.emplace_back(sample, weight);
            weights_sum += weight;
        }

        const std::vector<Sample>& list() const {
            return samples;
        }

        double sum() const {
            return weights_sum;
        }
    };
}

namespace {
    typedef std::mt19937 Random;

    using vcf::HOM;
    using vcf::HET;
    using vcf::HOMREF;
    using vcf::MISSING;
    using vcf::AlleleType;
    using vcf::Features;
    using vcf::Labels;
    using vcf::NodePtr;
    using vcf::Bags;

    double variance(const std::vector<double>& weights) {
        // Jeffrey's prior
        assert(weights.size() == 3);
        double sum_alpha = std::accumulate(weights.begin(), weights.end(), 0.0) + weights.size() / 2.0;
        std::vector<double> alpha(weights.size());
        std::transform(weights.begin(), weights.end(), alpha.begin(), [&sum_alpha](double x) {
            return (x + 0.5) / sum_alpha;
        });
        // linear model
        double mean = 0.0;
        for (size_t i = 0; i < weights.size(); i++) {
            mean += i * alpha[i];
        }
        double error = 0.0;
        for (size_t i = 0; i < weights.size(); i++) {
            error += alpha[i] * (i - mean) * (i - mean);
        }
        return error;
    }

    std::pair<double, double> weight_sums(const std::vector<double>& left_weights,
                                          const std::vector<double>& right_weights) {
        double left_sum = std::accumulate(left_weights.begin(), left_weights.end(), 0.0);
        double right_sum = std::accumulate(right_weights.begin(), right_weights.end(), 0.0);
        return {left_sum, right_sum};
    }

    class InnerNode : public vcf::Node {
        NodePtr left_node;
        NodePtr right_node;
        int var;
        vcf::AlleleType separator;

    public:
        InnerNode(std::vector<double>&& class_weights, NodePtr& left_node, NodePtr& right_node, vcf::AlleleType sep,
                  int variable);
        double predict(std::vector<vcf::AlleleType>& features) const override;
        static double joint_accuracy(NodePtr& left_node, NodePtr& right_node);
    };

    class LeafNode : public vcf::Node {
    public:
        explicit LeafNode(std::vector<double>&& class_weights);
        double predict(std::vector<vcf::AlleleType>& features) const override;
    };


    InnerNode::InnerNode(std::vector<double>&& class_weights, NodePtr& left_node, NodePtr& right_node,
                         vcf::AlleleType sep, int variable)
            :Node(std::move(class_weights)), left_node(left_node), right_node(right_node), var(variable),
             separator(sep) {
        acc = joint_accuracy(left_node, right_node);
    }

    double InnerNode::predict(std::vector<vcf::AlleleType>& features) const {
        vcf::AlleleType allele = features[var];
        if (allele != vcf::MISSING) {
            if (to_int(allele) <= to_int(separator)) {
                return left_node->predict(features);
            } else {
                return right_node->predict(features);
            }
        } else {
            double sum = std::accumulate(class_weights.begin(), class_weights.end(), 0.0);
            double left_ratio = class_weights[to_int(vcf::HOMREF)] / sum;
            double right_ratio = class_weights[to_int(vcf::HOM)] / sum;
            double het_ratio = class_weights[to_int(vcf::HET)] / sum;
            if (separator == vcf::HET) {
                left_ratio += het_ratio;
            } else {
                right_ratio += het_ratio;
            }
            return left_ratio * left_node->predict(features) + right_ratio * right_node->predict(features);
        }
    }

    double InnerNode::joint_accuracy(NodePtr& left_node, NodePtr& right_node) {
        auto sums = weight_sums(left_node->weights(), right_node->weights());
        double sum = sums.first + sums.second;
        return  (sums.first / sum) * left_node->accuracy() + (sums.second / sum) * right_node->accuracy();
    }

    LeafNode::LeafNode(std::vector<double>&& class_weights) :Node(std::move(class_weights)) {
        acc = variance(weights());
    }

    double LeafNode::predict(std::vector<vcf::AlleleType>& features) const {
        return prediction(class_weights);
    }

    class Counts {
        double counts[3] = {0, 0, 0};
    public:
        Counts() = default;

        Counts(double ref, double het, double alt) :counts{ref, het, alt} {}

        void add(AlleleType type, double weight) {
            counts[type] += weight;
        }

        double ref() const {
            return counts[0];
        }

        double het() const {
            return counts[1];
        }

        double alt() const {
            return counts[2];
        }

        double hom_ratio() const {
            return ref() / sum();
        }

        double het_ratio() const {
            return het() / sum();
        }

        double alt_ratio() const {
            return alt() / sum();
        }

        double entropy() const {
            double ratios[3] = {hom_ratio(), het_ratio(), alt_ratio()};
            double ret = 0.0;
            for (double r: ratios) {
                if (r > vcf::DecisionTree::EPS) {
                    ret -= r * std::log(r);
                }
            }
            return ret;
        }

        double sum() const {
            return ref() + het() + alt();
        }
    };

    Counts operator+(const Counts& one, const Counts& another) {
        return {one.ref() + another.ref(), one.het() + another.het(), one.alt() + another.alt()};
    }

    class Split {
        Bags l;
        Bags r;
        double gain;
    public:
        Split(Bags&& left, Bags&& right, double gain) :l(std::move(left)), r(std::move(right)), gain(gain) {}
        Split(Split&& other) :l(std::move(other.l)), r(std::move(other.r)), gain(other.gain) {}
        const Bags& left() {
            return l;
        }
        const Bags& right() {
            return r;
        }

        double score() {
            return gain;
        }
    };

    Counts counts(const Bags& bags, const Labels& labels) {
        Counts ret{};
        for (auto s: bags.list()) {
            ret.add(labels[s.sample()], s.weight());
        }
        return ret;
    }

    Split split(const Bags& curr, AlleleType splitBy, const std::vector<AlleleType>& features,
                                const Labels& labels, bool mock = true) {
        Bags left, right;
        Counts left_nm, right_nm;
        Counts cnts = counts(curr, labels);

        auto list = curr.list();
        double left_ratio = cnts.hom_ratio();
        if (splitBy == HET) {
            left_ratio += cnts.het_ratio();
        }

        for (auto& el: list) {
            auto allele = features.at(el.sample());
            if (allele == MISSING) {
                if (!mock) {
                    left.add(el.sample(), el.weight() * left_ratio);
                    right.add(el.sample(), el.weight() * (1.0 - left_ratio));
                }
            } else {
                if (allele <= splitBy) {
                    if (!mock) {
                        left.add(el.sample(), el.weight());
                    }
                    left_nm.add(labels[el.sample()], el.weight());
                } else {
                    if (!mock) {
                        right.add(el.sample(), el.weight());
                    }
                    right_nm.add(labels[el.sample()], el.weight());
                }
            }
        }
        Counts all_nm = left_nm + right_nm;
        double nm_ratio = all_nm.sum() / cnts.sum();
        double left_ratio_nm = left_nm.sum() / all_nm.sum();
        double split_entropy = left_ratio_nm * left_nm.entropy() + (1.0 - left_ratio_nm) * right_nm.entropy();
        double gain = nm_ratio * (all_nm.entropy() - split_entropy) - vcf::DecisionTree::EPS;
        return {std::move(left), std::move(right), gain};
    }

    std::vector<int> sample(size_t n, size_t k, Random& random) {
        assert(k <= n);
        std::vector<int> arr(n);
        std::iota(arr.begin(), arr.end(), 0);
        for (size_t i = 0; i < k; i++) {
            std::uniform_int_distribution<int> positions(0, arr.size() - i - 1);
            std::swap(arr[i], arr[positions(random)]);
        }
        return {arr.begin(), arr.begin() + k};
    }

    NodePtr prune(NodePtr left, NodePtr right, std::vector<double>&& class_weights, AlleleType sep,
                                int variable) {
        double common_variance = variance(class_weights);
        double joint_variance = InnerNode::joint_accuracy(left, right);
        if (common_variance < joint_variance - vcf::DecisionTree::EPS) {
            return std::make_shared<LeafNode>(std::move(class_weights));
        } else {
            return std::make_shared<InnerNode>(std::move(class_weights), left, right, sep, variable);
        }
    }
}

namespace vcf {
    double DecisionTree::predict(std::vector<AlleleType>& features) const {
        Bags bags;
        bags.add(0, 1.0);
        double ret = root->predict(features);
        if (ret < -EPS || ret > 2.0 + EPS) {
            throw std::logic_error("Error: predicted genotype is out of range");
        }
        return ret;
    }

    DecisionTree::DecisionTree(NodePtr root) :root(std::move(root)) {}

    DecisionTree::DecisionTree(DecisionTree&& other) noexcept :root(std::move(other.root)){}

    double DecisionTree::accuracy() {
        return root->accuracy();
    }

    Node::Node(std::vector<double>&& class_weights) :class_weights(std::move(class_weights)){}

    std::vector<double> Node::weights() {
        return class_weights;
    }

    double Node::prediction(const std::vector<double>& alpha) {
        // Dir(1,1,1) prior
        double sum = std::accumulate(alpha.begin(), alpha.end(), 0.0);
        if (sum < EPS) {
            sum = EPS;
        }

        std::vector<double> rel_alpha;
        std::for_each(alpha.begin(), alpha.end(), [&sum, &rel_alpha](double x){
            rel_alpha.push_back(x / sum);
        });
        return rel_alpha[1] + 2 * rel_alpha[2];
    }

    double Node::accuracy() {
        return acc;
    }

    TreeBuilder::TreeBuilder(const Features& features, const Labels& labels, size_t max_features)
        :features(features), values(labels), max_features(max_features) {}

    DecisionTree TreeBuilder::build_a_tree(Random& random, bool bagging) const {
        Bags tmp;
        for (size_t i = 0; i < values.size(); i++) {
            if (values[i] != MISSING) {
                tmp.add(i, 1.0);
            }
        }
        Bags bags;
        for (size_t i = 0; i < values.size(); i++) {
                switch(values[i]) {
                    case HOMREF: case HET: case HOM:
                        bags.add(i, 1.0);
                        break;
                    default:
                        continue;
                }
        }

        if (features.empty()) {
            auto cts = counts(tmp, values);
            std::vector<double> weights{cts.ref(), cts.het(), cts.alt()};
            return DecisionTree(std::make_shared<LeafNode>(std::move(weights)));
        }
        if (bagging) {
            Bags randomBags(bags, random);
            return DecisionTree(buildSubtree(randomBags, random));
        }
        return DecisionTree(buildSubtree(bags, random));
    }

    NodePtr TreeBuilder::buildSubtree(const Bags& bags, Random& random) const {
        auto vars = sample(features.size(), max_features, random);
        int var_best = -1;
        AlleleType best_split = MISSING;
        double best_score = 0.0;

        for (int var: vars) {
            auto hom_split = split(bags, HOMREF, features[var], values);
            auto het_split = split(bags, HET, features[var], values);
            if (hom_split.score() > best_score) {
                var_best = var;
                best_split = HOMREF;
                best_score = hom_split.score();
            }
            if (het_split.score() > best_score) {
                var_best = var;
                best_split = HET;
                best_score = het_split.score();
            }
        }

        auto cnts = counts(bags, values);
        std::vector<double> cs{cnts.ref(), cnts.het(), cnts.alt()};
        if (best_split == MISSING) {
            return std::make_shared<LeafNode>(std::move(cs));
        } else {
            auto the_best_split_ever = split(bags, best_split, features[var_best], values, false);
            auto& left = the_best_split_ever.left();
            auto& right = the_best_split_ever.right();
            auto left_subtree = buildSubtree(left, random);
            auto right_subtree = buildSubtree(right, random);
            return prune(left_subtree, right_subtree, std::move(cs), best_split, var_best);
        }
    }

    RandomForest::RandomForest(const TreeBuilder& treeBuilder, cxxpool::thread_pool& pool, size_t ntrees) {
        std::vector<std::future<DecisionTree>> futures;
        for (size_t i = 0; i < ntrees; i++) {
            int seed = rand();
            futures.push_back(pool.push([seed, &treeBuilder]() -> DecisionTree {
                Random random(seed);
                return treeBuilder.build_a_tree(random);
            }));
        }
        for (size_t i = 0; i < ntrees; i++) {
            futures[i].wait();
        }
        for (size_t i = 0; i < ntrees; i++) {
            predictors.push_back(futures[i].get());
        }
    }

    double RandomForest::predict(std::vector<AlleleType>& features) {
        double sum = 0.0;
        for (DecisionTree& tree: predictors) {
            sum += tree.predict(features);
        }
        return sum / predictors.size();
    }
}
