#include "include/vcf_parser.h"
#include <Rcpp.h>
#include <boost/algorithm/string/predicate.hpp>
#include <fstream>
#include <iostream>

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
               const CharacterVector& binary_prefix, const NumericVector& missingRateThreshold) {
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
            gmatrix_handler.reset(new RGenotypeMatrixHandler(ss, vs, stats, missingRateThreshold[0]));
            parser.register_handler(gmatrix_handler, 1);
            if (predictMissing[0]) {
                predicting_handler = make_shared<PredictingHandler>(ss, *gmatrix_handler, 250000, 100);
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

namespace {
    struct Counts {
        vector<int> hom;
        vector<int> het;
        vector<int> alt;
    };

    int ceiling_devision(size_t x, size_t y) {
        return (x + y - 1) / y;
    }

    class CountsReader {
        const vector<size_t>& samples;
        const MemoryMappedScanner& scanner;
        int DP;
        int GQ;
        int total_samples;
    public:
        CountsReader(const vector<size_t>& samples, const MemoryMappedScanner& scanner, int DP, int GQ, int tot_samples)
            :samples(samples), scanner(scanner), DP(DP), GQ(GQ), total_samples(tot_samples) {}

        std::tuple<int, int, int> read(size_t position) const {
            int homref = 0, het = 0, hom = 0;
            for (size_t sample_position: samples) {
                Allele allele = BinaryAllele::toAllele(scanner.scan(position * total_samples + sample_position));
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
            return {homref, het, hom};
        }
    };

    Counts parallel_read(const vector<size_t> &positions, const CountsReader& reader) {
        const int threads = 16;

        cxxpool::thread_pool pool{threads};
        std::vector<std::future<Counts>> futures;
        int jobs_per_thread = ceiling_devision(positions.size(), threads);

        int curr_jobs = 0;
        vector<size_t> thread_jobs;
        for (int i = 0; i < positions.size(); i++) {
            thread_jobs.push_back(positions[i]);
            if (curr_jobs == jobs_per_thread || i == positions.size() - 1) {
                futures.push_back(pool.push([thread_jobs, &reader]() -> Counts {
                    Counts counts = {};
                    for (size_t pos: thread_jobs) {
                        auto cts = reader.read(pos);
                        counts.hom.push_back(std::get<0>(cts));
                        counts.het.push_back(std::get<1>(cts));
                        counts.alt.push_back(std::get<2>(cts));
                    }
                    return counts;
                }));
                curr_jobs = 0;
                thread_jobs = vector<size_t>();
            }
        }

        Counts ret = {};

        for (auto& future: futures) {
            Counts done = future.get();
            ret.hom.insert(ret.hom.end(), done.hom.begin(), done.hom.end());
            ret.het.insert(ret.het.end(), done.het.begin(), done.het.end());
            ret.alt.insert(ret.alt.end(), done.alt.begin(), done.alt.end());
        }
        return ret;
    }
}

// [[Rcpp::export]]
List parse_binary_file(const CharacterVector& variants, const CharacterVector& samples,
        const CharacterVector& binary_file, const CharacterVector& metafile,
        IntegerVector& requiredDP, IntegerVector requiredGQ) {
    try {
        std::cerr << std::endl << sizeof(vcf::BinaryAllele) << std::endl;
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
        size_t n = 0;
        for (; iss >> sample; n++) {
            sample_positions.insert({sample, n});
        }
        vector<size_t> positions;
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

        Rcpp::LogicalVector present_variants;
        vector<size_t> variant_pos;
        for (const char* var: variants) {
            auto variant = Variant::parseVariants(string(var))[0];
            auto it = variant_positions.find(variant);
            present_variants.push_back(it != variant_positions.end());
            if (it != variant_positions.end()) {
               variant_pos.push_back(it->second);
            }
        }

        CountsReader reader(positions, scanner, DP, GQ, n);
        Counts counts = parallel_read(variant_pos, reader);

        List ret;
        ret["variant"] = variants[present_variants];
        ret["HOM_REF"] = NumericVector(counts.hom.begin(), counts.hom.end());
        ret["HET"] = NumericVector(counts.het.begin(), counts.het.end());
        ret["HOM_ALT"] = NumericVector(counts.alt.begin(), counts.alt.end());
        ret["total"] = positions.size();

        return ret;
    } catch (ParserException& e) {
        ::Rf_error(e.get_message().c_str());
    }
}
