#include "include/vcf_parser.h"
#include <Rcpp.h>
#include <boost/algorithm/string/predicate.hpp>
#include <fstream>

#include "include/vcf_binary.h"
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
            Rf_warning(e.get_message().c_str());
        }
    public:
        using VCFParser::VCFParser;
    };

    class RGenotypeMatrixHandler: public GenotypeMatrixHandler {
    public:
        using GenotypeMatrixHandler::GenotypeMatrixHandler;

        NumericMatrix result() {
            NumericMatrix res(gmatrix.size(), samples.size());
            for (size_t i = 0; i < gmatrix.size(); i++) {
                for (size_t j = 0; j < samples.size(); j++) {
                    float val = gmatrix[i][j];
                    if (val == to_int(vcf::MISSING)) {
                        val = NA_REAL;
                    }
                    res[j * gmatrix.size() + i] = val;
                }
            }
            vector<string> row_names;
            for_each(variants.begin(), variants.end(), [&row_names](Variant& v){
                row_names.push_back((string)v);
            });
            rownames(res) = CharacterVector(row_names.begin(), row_names.end());
            return res;
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

vector<size_t> predictor_positions(const vector<string>& samples, const CharacterVector& exclude) {
    unordered_set<string> exclude_set;
    for_each(exclude.begin(), exclude.end(), [&exclude_set](const char* s) {
        string sample(s);
        exclude_set.insert(sample);
    });
    vector<size_t> ret;
    for (size_t i = 0; i < samples.size(); i++) {
        if (exclude_set.find(samples[i]) == exclude_set.end()) {
            ret.push_back(i);
        }
    }
    return ret;
}

// [[Rcpp::export]]
List parse_vcf(const CharacterVector& filename, const CharacterVector& samples,
               const CharacterVector& bad_positions, const CharacterVector& variants,
               const IntegerVector& DP, const IntegerVector& GQ, const LogicalVector& gmatrix,
               const LogicalVector& predictMissing, const CharacterVector& excludedPredictors,
               const CharacterVector& regions, const CharacterVector& binary_prefix) {
    List ret;
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
            gmatrix_handler.reset(new RGenotypeMatrixHandler(ss, vs, stats));
            parser.register_handler(gmatrix_handler, 1);
            if (predictMissing[0]) {
                auto predictor_samples = predictor_positions(ss, excludedPredictors);
                predicting_handler = make_shared<PredictingHandler>(ss, *gmatrix_handler, 250000, 100,
                                                                    predictor_samples);
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
        ::Rf_error(e.get_message().c_str());
    }
    return ret;
}

// [[Rcpp::export]]
List parse_binary_file(const CharacterVector& variants, const CharacterVector& samples,
        const CharacterVector& binary_file, const CharacterVector& metafile,
        IntegerVector& requiredDP, IntegerVector requiredGQ) {
    try {
        int DP = requiredDP[0];
        int GQ = requiredGQ[0];

        MemoryMappedScanner scanner((string) binary_file[0]);
        ifstream fin(metafile[0]);
        string line;
        getline(fin, line);
        istringstream iss(line);
        unordered_map<string, int> sample_positions;
        unordered_map<Variant, int> variant_positions;
        string sample;
        int n = 0;
        for (; iss >> sample; n++) {
            sample_positions.insert({sample, n});
        }
        vector<int> positions;
        for (const char *s: samples) {
            positions.push_back(sample_positions[string(s)]);
        }

        int i = 0;
        while (getline(fin, line)) {
            auto vars = Variant::parseVariants(line);
            for (Variant v: vars) {
                variant_positions.insert({v, i++});
            }
        }

        List ret;
        vector<int> homrefs, hets, homs;
        for (const char *var: variants) {
            int homref = 0, het = 0, hom = 0;
            int pos = variant_positions.at(Variant::parseVariants(string(var))[0]);
            for (int sample_position: positions) {
                Allele allele = BinaryAllele::toAllele(scanner.scan(pos * n + sample_position));
                if (allele.DP() >= DP && allele.GQ() >= GQ) {
                    switch (allele.alleleType()) {
                        case HOM:
                            ++hom;
                            break;
                        case HET:
                            ++het;
                            break;
                        case HOMREF:
                            ++homref;
                            break;
                        default:
                            break;
                    }
                }
            }
            homrefs.push_back(homref);
            hets.push_back(het);
            homs.push_back(hom);
        }
        ret["variant"] = variants;
        ret["HOM_REF"] = NumericVector(homrefs.begin(), homrefs.end());
        ret["HET"] = NumericVector(hets.begin(), hets.end());
        ret["HOM_ALT"] = NumericVector(homs.begin(), homs.end());
        ret["total"] = positions.size();
        return ret;
    } catch (ParserException& e) {
        ::Rf_error(e.get_message().c_str());
    }
}
