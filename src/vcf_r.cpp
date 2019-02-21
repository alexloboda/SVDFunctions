#include "vcf_parser.h"
#include <Rcpp.h>
#include <boost/algorithm/string/predicate.hpp>
#include <iostream>
#include <fstream>
#include "zstr/zstr.hpp"
#include "zstr/strict_fstream.hpp"

namespace {
    using namespace Rcpp;
    using namespace vcf;
    using namespace std;
    using boost::algorithm::ends_with;

    class Parser: public VCFParser {
        void handle_error(const vcf::ParserException& e) override {
            Rf_warning(e.get_message().c_str());
        }
    public:
        using VCFParser::VCFParser;
    };

    class RGenotypeMatrixHandler: public GenotypeMatrixHandler {
    public:
        using GenotypeMatrixHandler::GenotypeMatrixHandler;

        IntegerMatrix result() {
            IntegerMatrix res(gmatrix.size(), samples.size());
            for (int i = 0; i < gmatrix.size(); i++) {
                for (int j = 0; j < samples.size(); j++) {
                    int val = gmatrix[i][j];
                    if (val == vcf::MISSING) {
                        val = NA_INTEGER;
                    }
                    res[i * samples.size() + j] = val;
                }
            }
            return res;
        }
    };

    class RCallRateHandler: public CallRateHandler {
    public:
        using CallRateHandler::CallRateHandler;

        NumericMatrix result() {
            NumericMatrix result(ranges.size(), samples.size());
            for (int i = 0; i < ranges.size(); i++) {
                for (int j = 0; j < samples.size(); j++) {
                    result[i * samples.size() + j] = (double)call_rate_matrix[i][j] / variants;
                }
            }
            return result;
        }
    };
}

VCFFilter filter(const CharacterVector& samples, const CharacterVector& bad_positions,
        const CharacterVector& variants, int DP, int GQ) {
    VCFFilter filter(DP, GQ);

    if (samples.length() > 0) {
        vector<string> ss;
        transform(samples.begin(), samples.end(), ss.begin(), [](const char *s) { return string(s); });
        filter.add_samples(ss);
    }

    if (bad_positions.length() > 0) {
        vector<Position> bads;
        transform(bad_positions.begin(), bad_positions.end(), bads.begin(), [](const char *s) {
            return Position::parse_position(string(s));
        });
        filter.add_bad_variants(bads);
    }

    if (variants.length() > 0) {
        vector<Variant> vs;
        for_each(variants.begin(), variants.end(), [&vs](const char *s) {
            vector<Variant> variants = Variant::parseVariants(string(s));
            vs.insert(vs.end(), variants.begin(), variants.end());
        });
        filter.set_available_variants(vs);
    }
    return filter;
}

vector<vcf::Range> parse_regions(const CharacterVector& regions){
    vector<vcf::Range> ranges;
    transform(regions.begin(), regions.end(), ranges.begin(), [](const char* str){
        return vcf::Range::parseRange(string(str));
    });
}

// [[Rcpp::export]]
List parse_vcf(const CharacterVector& filename, const CharacterVector& samples,
               const CharacterVector& bad_positions, const CharacterVector& allowed_variants,
               const IntegerVector& DP, const IntegerVector& GQ, const CharacterVector& regions,
               const LogicalVector& ret_gmatrix, const CharacterVector& binary_prefix) {
    const char* name = filename[0];
    unique_ptr<std::istream> in(new zstr::ifstream(name));

    Parser parser(*in, filter(samples, bad_positions, allowed_variants, DP[0], GQ[0]));
    parser.parse_header();
    auto ss = parser.sample_names();
    shared_ptr<RGenotypeMatrixHandler> gmatrix_handler;
    shared_ptr<BinaryFileHandler> binary_handler;
    shared_ptr<RCallRateHandler> callrate_handler;

    if (ret_gmatrix[0]) {
        gmatrix_handler.reset(new RGenotypeMatrixHandler(ss));
        parser.register_handler(gmatrix_handler);
    }

    if (regions.length() > 0) {
        callrate_handler.reset(new RCallRateHandler(ss, parse_regions(regions)));
        parser.register_handler(callrate_handler);
    }

    if (binary_prefix.length() > 0) {
        string prefix = string(binary_prefix[0]);
        binary_handler.reset(new BinaryFileHandler(ss, prefix + "_bin", prefix + "_meta"));
        parser.register_handler(binary_handler);
    }

    parser.parse_genotypes();
    List ret;
    vector<const char*> samples_cstr;
    transform(ss.begin(), ss.end(), samples_cstr.begin(), [](string& s){return s.c_str();});
    ret["samples"] = CharacterVector(samples_cstr.begin(), samples_cstr.end());
    if (ret_gmatrix[0]){
        ret["genotype"] = gmatrix_handler->result();
    }
    if (regions.length() > 0) {
        ret["callrate"] = callrate_handler->result();
    }
    return ret;
}
