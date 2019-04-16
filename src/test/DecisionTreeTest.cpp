#define CATCH_CONFIG_MAIN

#include "../include/test/catch.hpp"
#include "../include/genotype_predictor.h"

namespace {
    using namespace vcf;
    const double EPS = 1e-8;
}

TEST_CASE("Decision trees constructed correctly"){
    Features data{{HOM, HET, HOM, HET, HOM, HET},
                  {HOM, HOM, HOM, HET, HET, HET}};

    Labels labels{HOMREF, HOMREF, HET, HET, HOM, HOM};

    SECTION("manual accuracy test") {
        TreeBuilder tree_builder(data, labels, 2);
        std::mt19937 random{42};
        DecisionTree tree = tree_builder.build_a_tree(random, false);
        REQUIRE(fabs(tree.accuracy() - 5.0 / 9.0) < EPS);
    }

    SECTION("random forest check") {
        TreeBuilder tree_builder(data, labels, 2);
        RandomForest rf{tree_builder, 4};
        std::vector<AlleleType> features{HOM, HOM};
        REQUIRE(rf.predict(features) <= 1);
    }
}
