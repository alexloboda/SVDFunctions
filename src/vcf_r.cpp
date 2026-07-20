#include "include/vcf_parser.h"
#include <Rcpp.h>
#include <boost/algorithm/string/predicate.hpp>
#include <fstream>

#include "include/third-party/zstr/zstr.hpp"
#include "include/third-party/zstr/strict_fstream.hpp"
#include "include/vcf_stats.h"
#include "include/vcf_predicting_handler.h"

namespace {
    using namespace Rcpp;
    using namespace vcf;
    using namespace std;
    using boost::algorithm::ends_with;

    class Parser: public VCFParser {
        void handle_error(const vcf::ParserException& e) override {
            Rcpp::warning(e.get_message());
        }
    public:
        using VCFParser::VCFParser;
    };

    class RGenotypeMatrixHandler: public GenotypeMatrixHandler {
    public:
        using GenotypeMatrixHandler::GenotypeMatrixHandler;

        List result() {
            NumericMatrix res(gmatrix.size(), samples.size());
            LogicalMatrix predicted(missing.size(), samples.size());
            for (size_t i = 0; i < gmatrix.size(); i++) {
                for (size_t j = 0; j < samples.size(); j++) {
                    float val = gmatrix[i][j];
                    if (val == to_int(vcf::MISSING)) {
                        val = NA_REAL;
                    }
                    res[j * gmatrix.size() + i] = val;
                    predicted[j * missing.size() + i] = missing[i][j];
                }
            }
            vector<string> row_names;
            for_each(variants.begin(), variants.end(), [&row_names](Variant& v){
                row_names.push_back((string)v);
            });
            rownames(res) = CharacterVector(row_names.begin(), row_names.end());
            rownames(predicted) = CharacterVector(row_names.begin(), row_names.end());
            List ret;
            ret["genotype"] = res;
            ret["predicted"] = predicted;
            return ret;
        }
    };

    class RCallRateHandler: public CallRateHandler {
    public:
        using CallRateHandler::CallRateHandler;

        NumericMatrix result() {
            vector<string> non_empty;
            for (size_t i = 0; i < ranges.size(); i++) {
                if (n_variants[i] > 0) {
                    non_empty.push_back((std::string)ranges[i]);
                }
            }
            NumericMatrix result(non_empty.size(), samples.size());
            int curr = 0;
            for (size_t i = 0; i < ranges.size(); i++) {
                if (n_variants[i] == 0) {
                    continue;
                }
                for (size_t j = 0; j < samples.size(); j++) {
                    result[j * non_empty.size() + curr] = (double)call_rate_matrix[i][j] / n_variants[i];
                }
                ++curr;
            }
            rownames(result) = CharacterVector(non_empty.begin(), non_empty.end());
            return result;
        }
    };
}

VCFFilter filter(const CharacterVector& samples, const CharacterVector& bad_positions,
        int DP, int GQ) {
    VCFFilter filter(DP, GQ);

    if (samples.length() > 0) {
        vector<string> ss;
        for_each(samples.begin(), samples.end(), [&ss](const char *s) { ss.emplace_back(s); });
        filter.add_samples(ss);
    }

    if (bad_positions.length() > 0) {
        vector<Position> bads;
        for_each(bad_positions.begin(), bad_positions.end(), [&bads](const char *s) {
            bads.push_back(Position::parse_position(string(s)));
        });
        filter.add_bad_variants(bads);
    }

    return filter;
}

vector<vcf::Range> parse_regions(const CharacterVector& regions){
    vector<vcf::Range> ranges;
    for_each(regions.begin(), regions.end(), [&ranges](const char* str){
        ranges.push_back(vcf::Range::parseRange(string(str)));
    });
    return ranges;
}

// [[Rcpp::export]]
List parse_vcf(const CharacterVector& filename, const CharacterVector& samples,
               const CharacterVector& bad_positions, const CharacterVector& variants,
               const IntegerVector& DP, const IntegerVector& GQ, const LogicalVector& gmatrix,
               const LogicalVector& predictMissing, const CharacterVector& regions,
               const CharacterVector& binary_prefix, const NumericVector& missingRateThreshold,
               Rcpp::Nullable<int> seed) {
    List ret;
    unsigned int random_seed = Rcpp::as<int>(seed);
    try {
        const char *name = filename[0];
        unique_ptr<std::istream> in(new zstr::ifstream(name));

        VCFFilterStats stats;
        Parser parser(*in, filter(samples, bad_positions, DP[0], GQ[0]), stats);
        parser.parse_header();
        auto ss = parser.sample_names();
        shared_ptr<RGenotypeMatrixHandler> gmatrix_handler;
        shared_ptr<BinaryFileHandler> binary_handler;
        shared_ptr<RCallRateHandler> callrate_handler;
        shared_ptr<PredictingHandler> predicting_handler;

        if (gmatrix[0]) {
            vector<Variant> vs;
            for_each(variants.begin(), variants.end(), [&vs](const char *s) {
                vector<Variant> variants = Variant::parseVariants(string(s));
                vs.insert(vs.end(), variants.begin(), variants.end());
            });
            gmatrix_handler.reset(new RGenotypeMatrixHandler(ss, vs, stats, missingRateThreshold[0]));
            parser.register_handler(gmatrix_handler, 1);
            if (predictMissing[0]) {
                predicting_handler = make_shared<PredictingHandler>(ss, *gmatrix_handler, 250000, 100, random_seed);
                parser.register_handler(predicting_handler, 2);
            }
        }

        if (regions.length() > 0) {
            callrate_handler.reset(new RCallRateHandler(ss, parse_regions(regions)));
            parser.register_handler(callrate_handler, 1);
        }

        if (binary_prefix.length() > 0) {
            string prefix = string(binary_prefix[0]);
            binary_handler.reset(new BinaryFileHandler(ss, prefix + "_bin", prefix + "_meta"));
            parser.register_handler(binary_handler, 1);
        }

        if (gmatrix_handler != nullptr || binary_handler != nullptr || callrate_handler != nullptr) {
            parser.parse_genotypes();
        }
        ret["samples"] = CharacterVector(ss.begin(), ss.end());
        if (gmatrix[0]) {
            if (predictMissing[0]) {
                predicting_handler->cleanup();
            }
            ret["genotype"] = gmatrix_handler->result();
        }
        if (regions.length() > 0) {
            ret["callrate"] = callrate_handler->result();
        }
        List ret_stats;
        for (Stat stat: vcf::statsList()) {
            ret_stats[to_string(stat)] = stats.value(stat);
        }
        ret["stats"] = ret_stats;
    } catch (ParserException& e) {
        Rcpp::stop(e.get_message());
    }
    return ret;
}