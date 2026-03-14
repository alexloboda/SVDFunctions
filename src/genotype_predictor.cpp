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
            const auto& list = bags.list();
            prefix_weights.reserve(list.size());
            samples.reserve(list.size());
            for (const auto& s: list) {
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

        void reserve(size_t n) {
            samples.reserve(n);
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

    Counts operator-(const Counts& one, const Counts& another) {
        return {one.ref() - another.ref(), one.het() - another.het(), one.alt() - another.alt()};
    }

    class Split {
        Bags l;
        Bags r;
        Counts left_counts_value;
        Counts right_counts_value;
        double gain;
    public:
        Split(Bags&& left, Bags&& right, Counts&& left_counts, Counts&& right_counts, double gain)
            :l(std::move(left)), r(std::move(right)), left_counts_value(std::move(left_counts)),
             right_counts_value(std::move(right_counts)), gain(gain) {}
        Split(Split&& other)
            :l(std::move(other.l)), r(std::move(other.r)),
             left_counts_value(std::move(other.left_counts_value)),
             right_counts_value(std::move(other.right_counts_value)), gain(other.gain) {}
        const Bags& left() {
            return l;
        }
        const Bags& right() {
            return r;
        }
        const Counts& left_counts() const {
            return left_counts_value;
        }
        const Counts& right_counts() const {
            return right_counts_value;
        }

        double score() const {
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

    double split_gain(const Counts& total_counts, const Counts& left_counts, const Counts& right_counts) {
        Counts observed_counts = left_counts + right_counts;
        double nm_ratio = observed_counts.sum() / total_counts.sum();
        double left_ratio_nm = left_counts.sum() / observed_counts.sum();
        double split_entropy = left_ratio_nm * left_counts.entropy() + (1.0 - left_ratio_nm) * right_counts.entropy();
        double gain = nm_ratio * (observed_counts.entropy() - split_entropy) - vcf::DecisionTree::EPS;
        return gain;
    }

    struct SplitScores {
        double homref_gain;
        double het_gain;
    };

    SplitScores evaluate_split_scores(const Bags& curr, const Counts& total_counts,
                                      const std::vector<AlleleType>& feature_values,
                                      const Labels& labels) {
        Counts left_homref;
        Counts right_homref;
        Counts left_het;
        Counts right_het;

        for (const auto& el : curr.list()) {
            auto allele = feature_values[el.sample()];
            if (allele == MISSING) {
                continue;
            }

            auto label = labels[el.sample()];
            if (allele <= HOMREF) {
                left_homref.add(label, el.weight());
                left_het.add(label, el.weight());
            } else if (allele <= HET) {
                right_homref.add(label, el.weight());
                left_het.add(label, el.weight());
            } else {
                right_homref.add(label, el.weight());
                right_het.add(label, el.weight());
            }
        }

        return {
            split_gain(total_counts, left_homref, right_homref),
            split_gain(total_counts, left_het, right_het)
        };
    }

    Split split(const Bags& curr, const Counts& total_counts, AlleleType splitBy,
                const std::vector<AlleleType>& feature_values, const Labels& labels) {
        Bags left, right;
        Counts left_nm, right_nm;
        Counts left_total, right_total;
        const auto& list = curr.list();
        left.reserve(list.size());
        right.reserve(list.size());

        double left_ratio = total_counts.hom_ratio();
        if (splitBy == HET) {
            left_ratio += total_counts.het_ratio();
        }

        for (const auto& el: list) {
            auto allele = feature_values[el.sample()];
            if (allele == MISSING) {
                double left_weight = el.weight() * left_ratio;
                double right_weight = el.weight() * (1.0 - left_ratio);
                left.add(el.sample(), left_weight);
                right.add(el.sample(), right_weight);
                left_total.add(labels[el.sample()], left_weight);
                right_total.add(labels[el.sample()], right_weight);
            } else {
                if (allele <= splitBy) {
                    left.add(el.sample(), el.weight());
                    left_nm.add(labels[el.sample()], el.weight());
                    left_total.add(labels[el.sample()], el.weight());
                } else {
                    right.add(el.sample(), el.weight());
                    right_nm.add(labels[el.sample()], el.weight());
                    right_total.add(labels[el.sample()], el.weight());
                }
            }
        }
        double gain = split_gain(total_counts, left_nm, right_nm);
        return {std::move(left), std::move(right), std::move(left_total), std::move(right_total), gain};
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

    NodePtr build_subtree_impl(const Bags& bags, const Counts& total_counts, const Features& features,
                               const Labels& values, size_t max_features, Random& random) {
        auto vars = sample(features.size(), max_features, random);
        int var_best = -1;
        AlleleType best_split = MISSING;
        double best_score = 0.0;

        for (int var: vars) {
            auto scores = evaluate_split_scores(bags, total_counts, features[var], values);
            if (scores.homref_gain > best_score) {
                var_best = var;
                best_split = HOMREF;
                best_score = scores.homref_gain;
            }
            if (scores.het_gain > best_score) {
                var_best = var;
                best_split = HET;
                best_score = scores.het_gain;
            }
        }

        std::vector<double> cs{total_counts.ref(), total_counts.het(), total_counts.alt()};
        if (best_split == MISSING) {
            return std::make_shared<LeafNode>(std::move(cs));
        }

        auto best = split(bags, total_counts, best_split, features[var_best], values);
        auto& left = best.left();
        auto& right = best.right();
        auto left_subtree = build_subtree_impl(left, best.left_counts(), features, values, max_features, random);
        auto right_subtree = build_subtree_impl(right, best.right_counts(), features, values, max_features, random);
        return prune(left_subtree, right_subtree, std::move(cs), best_split, var_best);
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

    const std::vector<double>& Node::weights() const {
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

    double Node::accuracy() const {
        return acc;
    }

    TreeBuilder::TreeBuilder(const Features& features, const Labels& labels, size_t max_features)
        :features(features), values(labels), max_features(max_features) {}

    std::pair<DecisionTree, std::vector<unsigned char>> TreeBuilder::build_a_tree_with_inbag(Random& random, bool bagging) const {
        Bags tmp;
        Bags bags;
        for (size_t i = 0; i < values.size(); i++) {
            auto value = values[i];
            if (value != MISSING) {
                tmp.add((int)i, 1.0);
            }
            switch(value) {
                case HOMREF: case HET: case HOM:
                    bags.add((int)i, 1.0);
                    break;
                default:
                    continue;
            }
        }

        std::vector<unsigned char> inbag(values.size(), 0);

        if (features.empty()) {
            auto cts = counts(tmp, values);
            std::vector<double> weights{cts.ref(), cts.het(), cts.alt()};
            return {DecisionTree(std::make_shared<LeafNode>(std::move(weights))), std::move(inbag)};
        }

        if (bagging) {
            Bags randomBags(bags, random);
            for (const auto& s : randomBags.list()) {
                if (s.sample() >= 0 && (size_t)s.sample() < inbag.size()) {
                    inbag[(size_t)s.sample()] = 1;
                }
            }
            return {DecisionTree(buildSubtree(randomBags, random)), std::move(inbag)};
        }

        // No bagging: mark all eligible samples as in-bag.
        for (const auto& s : bags.list()) {
            if (s.sample() >= 0 && (size_t)s.sample() < inbag.size()) {
                inbag[(size_t)s.sample()] = 1;
            }
        }
        return {DecisionTree(buildSubtree(bags, random)), std::move(inbag)};
    }

    DecisionTree TreeBuilder::build_a_tree(Random& random, bool bagging) const {
        return build_a_tree_with_inbag(random, bagging).first;
    }

    NodePtr TreeBuilder::buildSubtree(const Bags& bags, Random& random) const {
        auto total_counts = counts(bags, values);
        return build_subtree_impl(bags, total_counts, features, values, max_features, random);
    }

    RandomForest::RandomForest(const TreeBuilder& treeBuilder, cxxpool::thread_pool& pool, size_t ntrees, unsigned int seed) {
        std::vector<std::future<std::pair<DecisionTree, std::vector<unsigned char>>>> futures;
        futures.reserve(ntrees);
        predictors.reserve(ntrees);
        inbag_masks.reserve(ntrees);
        for (size_t i = 0; i < ntrees; i++) {
            int tree_seed = seed + i;
            futures.push_back(pool.push([tree_seed, &treeBuilder]() -> std::pair<DecisionTree, std::vector<unsigned char>> {
                Random random(tree_seed);
                return treeBuilder.build_a_tree_with_inbag(random);
            }));
        }
        for (size_t i = 0; i < ntrees; i++) {
            futures[i].wait();
        }
        for (size_t i = 0; i < ntrees; i++) {
            auto built = futures[i].get();
            predictors.push_back(std::move(built.first));
            inbag_masks.push_back(std::move(built.second));
        }
    }

    double RandomForest::predict(std::vector<AlleleType>& features) {
        double sum = 0.0;
        for (DecisionTree& tree: predictors) {
            sum += tree.predict(features);
        }
        return sum / predictors.size();
    }

    double RandomForest::predict_oob(std::vector<AlleleType>& features, std::size_t sample_index) {
        if (predictors.empty()) {
            throw std::logic_error("Error: random forest has no predictors");
        }
        double sum = 0.0;
        std::size_t cnt = 0;
        for (std::size_t t = 0; t < predictors.size(); t++) {
            if (t >= inbag_masks.size()) {
                break;
            }
            if (sample_index >= inbag_masks[t].size()) {
                break;
            }
            if (inbag_masks[t][sample_index] == 0) {
                sum += predictors[t].predict(features);
                ++cnt;
            }
        }
        if (cnt == 0) {
            return predict(features);
        }
        return sum / (double)cnt;
    }
}
