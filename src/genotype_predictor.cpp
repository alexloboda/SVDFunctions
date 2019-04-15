#include <utility>

#include "include/genotype_predictor.h"
#include "include/vcf_primitives.h"

namespace {
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
}

namespace vcf {

    class Bags {
        std::vector<int> bags;
        std::vector<double> weights;
        double weights_sum = 0.0;

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

    double variance(const std::vector<double>& weights) {
        // uniform prior (beta(1, 1, 1))
        std::vector<double> alpha(weights.size());
        std::transform(weights.begin(), weights.end(), alpha.begin(), [](double x) { return x + 1; });
        double sum_alpha = std::accumulate(alpha.begin(), alpha.end(), 0.0);
        std::vector<double> rel_alpha;
        std::vector<double> variance;
        std::for_each(alpha.begin(), alpha.end(), [&](double x) {
            double alpha_hat = x / sum_alpha;
            rel_alpha.push_back(alpha_hat);
            variance.push_back((alpha_hat * (1 - alpha_hat)) / (sum_alpha + 1));
        });

        double cov_bc = (-rel_alpha[1] * rel_alpha[2]) / (sum_alpha + 1);
        return variance[1] + 4 * variance[2] + 4 * cov_bc;
    }

    std::pair<double, double> weight_sums(const std::vector<double>& left_weights,
                                          const std::vector<double>& right_weights) {
        double left_sum = std::accumulate(left_weights.begin(), left_weights.end(), 0.0);
        double right_sum = std::accumulate(right_weights.begin(), right_weights.end(), 0.0);
        return {left_sum, right_sum};
    }

    double joint_variance(std::vector<double>& left_weights, std::vector<double>& right_weights) {
        auto sums = weight_sums(left_weights, right_weights);
        double sep_variance = sums.first * variance(left_weights) + sums.second * variance(right_weights);
        return sep_variance / sums.first + sums.second;
    }

    class InnerNode : public vcf::Node {
        NodePtr left_node;
        NodePtr right_node;
        int var;
        vcf::AlleleType separator;

    public:
        InnerNode(std::vector<double>&& class_weights, NodePtr& left_node, NodePtr& right_node, vcf::AlleleType sep,
                  int variable);
        double predict(std::vector<vcf::AlleleType>& features) override;
        double accuracy() const override;
    };

    class LeafNode : public vcf::Node {
    public:
        explicit LeafNode(std::vector<double>&& class_weights);
        double predict(std::vector<vcf::AlleleType>& features) override;
        double accuracy() const override;
    };

    InnerNode::InnerNode(std::vector<double>&& class_weights, NodePtr& left_node, NodePtr& right_node,
                         vcf::AlleleType sep, int variable)
            :Node(std::move(class_weights)), left_node(left_node), right_node(right_node), var(variable),
             separator(sep) {}

    double InnerNode::predict(std::vector<vcf::AlleleType>& features) {
        vcf::AlleleType allele = features[var];
        if (allele != vcf::MISSING) {
            if (to_int(allele) <= to_int(separator)) {
                return left_node->predict(features);
            } else {
                return right_node->predict(features);
            }
        } else {
            double left_ratio = class_weights[to_int(vcf::HOMREF)];
            double right_ratio = class_weights[to_int(vcf::HOM)];
            double het_ratio = class_weights[to_int(vcf::HET)];
            if (separator == vcf::HET) {
                left_ratio += het_ratio;
            } else {
                right_ratio += het_ratio;
            }
            return left_ratio * left_node->predict(features) + right_ratio * right_node->predict(features);
        }
    }

    double InnerNode::accuracy() const {
        auto sums = weight_sums(left_node->weights(), right_node->weights());
        double sum = sums.first + sums.second;
        return (sums.first / sum) * left_node->accuracy() + (sums.second / sum) * right_node->accuracy();
    }

    LeafNode::LeafNode(std::vector<double>&& class_weights) :Node(std::move(class_weights)) {}

    double LeafNode::predict(std::vector<vcf::AlleleType>& features) {
        return prediction(class_weights);
    }

    double LeafNode::accuracy() const {
        return variance(class_weights);
    }

    typedef std::pair<Bags, Bags> Split;

    std::tuple<double, double, double> counts(Bags& bags, const Labels& labels) {
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
                case MISSING:
                    throw std::logic_error("Predictable values must not be NAs");
            }
        }
        return {hom, het, alt};
    }

    std::tuple<double, double, double> ratios(Bags& bags, const Labels& labels) {
        auto count = counts(bags, labels);

        double hom = std::get<0>(count);
        double het = std::get<1>(count);
        double alt = std::get<2>(count);

        double sum = hom + het + alt;
        return {hom / sum, het / sum, alt / sum};
    }

    std::pair<Bags, Bags> split(Bags& curr, AlleleType splitBy, const std::vector<AlleleType>& features,
                                const Labels& labels) {
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

    double score(Bags& samples, const Labels& labels) {
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

    double split_score(Split& split, const Labels& labels) {
        double sum = split.first.sum() + split.second.sum();
        double lr = split.first.sum() / sum;
        double rr = split.second.sum() / sum;
        return lr * score(split.first, labels) + rr * score(split.second, labels);
    }

    NodePtr prune(NodePtr left, NodePtr right, std::vector<double>&& class_weights, AlleleType sep,
                                int variable) {
        auto left_weights = left->weights();
        auto right_weights = right->weights();
        double sep_variance = joint_variance(left_weights, right_weights);
        double common_variance = variance(class_weights);
        if (common_variance < sep_variance - vcf::DecisionTree::EPS) {
            return std::make_shared<LeafNode>(std::move(class_weights));
        } else {
            return std::make_shared<InnerNode>(std::move(class_weights), left, right, sep, variable);
        }
    }
}

namespace vcf {
    double DecisionTree::predict(std::vector<AlleleType>& features) {
        Bags bags;
        bags.add(0, 1.0);
        return root->predict(features);
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

    double Node::prediction(std::vector<double>& alpha) {
        // Beta(1,1,1) prior
        double sum = std::accumulate(alpha.begin(), alpha.end(), 0.0) + alpha.size();
        std::transform(alpha.begin(), alpha.end(), alpha.begin(), [&sum](double x){
            return x / sum;
        });
        return alpha[1] + 2 * alpha[2];
    }
    TreeBuilder::TreeBuilder(Features& features, Labels& labels, size_t max_features) :features(features), values(labels),
                                                                                       max_features(max_features){}

    DecisionTree TreeBuilder::build_a_tree(Random& random) {
        Bags bags(values.size(), random);
        return DecisionTree(buildSubtree(bags, random));
    }

    NodePtr TreeBuilder::buildSubtree(Bags& bags, Random& random) {
        auto vars = sample(features.size(), max_features, random);
        int var_best = -1;
        AlleleType best_split = MISSING;
        double best_score = score(bags, values) - vcf::DecisionTree::EPS;

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

        auto cnts = counts(bags, values);
        std::vector<double> cs{std::get<0>(cnts), std::get<1>(cnts), std::get<2>(cnts)};
        if (best_split == MISSING) {
            return std::make_shared<LeafNode>(std::move(cs));
        } else {
            auto the_best_split_ever = split(bags, best_split, features[var_best], values);
            auto left = the_best_split_ever.first;
            auto right = the_best_split_ever.second;
            auto left_subtree = buildSubtree(left, random);
            auto right_subtree = buildSubtree(right, random);
            return prune(left_subtree, right_subtree, std::move(cs), best_split, var_best);
        }
    }
}
