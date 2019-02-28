#include <Rcpp.h>
#include <fstream>
#include <unordered_set>
#include <string>

#include "include/vcf_primitives.h"

namespace {
    using Rcpp::CharacterVector;
    using Rcpp::IntegerVector;

    using std::string;
    using std::vector;
    using std::unordered_set;

    using vcf::Variant;
    using vcf::AlleleBinary;

    vector<string> parse_samples(std::istream& in) {
        vector<string> samples;
        string line;
        getline(in, line);
        std::istringstream iss(line);
        string sample;
        for (int i = 0; iss >> sample; i++) {
            samples.push_back(sample);
        }
        return samples;
    }

    vector<Variant> parse_variants(std::istream& in) {
        vector<Variant> variants;
        string line;
        while (getline(in, line)) {
            if (std::all_of(line.begin(), line.end(), isspace)) {
                continue;
            }
            auto vars = Variant::parseVariants(line);
            for (const Variant& v: vars) {
                variants.push_back(v);
            }
        }
        return variants;
    }

}
