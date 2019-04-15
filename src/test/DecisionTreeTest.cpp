#define CATCH_CONFIG_MAIN
#include "../include/test/catch.hpp"
#include "../include/genotype_predictor.h"

namespace {
    using namespace vcf;
    const double EPS = 1e-8;
}

TEST_CASE("Decision trees constructed correctly"){
    SECTION("manual accuracy test") {
        std::vector<std::vector<vcf::AlleleType>> data{{HOM, HET, HOM, HET, HOM, HET},
                                                       {HOM, HOM, HOM, HET, HET, HET}};

        std::vector<vcf::AlleleType> labels{HOMREF, HOMREF, HET, HET, HOM, HOM};
        TreeBuilder tree_builder(data, labels, 2);
        std::mt19937 random{42};
        DecisionTree tree = tree_builder.build_a_tree(random);
        REQUIRE(fabs(tree.accuracy() - 0.06746032) < EPS);
    }
}
