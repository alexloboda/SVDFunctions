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

namespace {
    class Counts {
        vector<int> hom;
        vector<int> het;
        vector<int> alt;
    public:
        void push(int homc, int hetc, int altc) {
            hom.push_back(homc);
            het.push_back(hetc);
            alt.push_back(altc);
        }

        void add(int entry, int homc, int hetc, int altc) {
            hom[entry] += homc;
            het[entry] += hetc;
            alt[entry] += altc;
        }

        void append(const Counts& counts) {
            hom.insert(hom.end(), counts.hom.begin(), counts.hom.end());
            het.insert(het.end(), counts.het.begin(), counts.het.end());
            alt.insert(alt.end(), counts.alt.begin(), counts.alt.end());
        }

        List get_table() {
            List ret;
            ret["hom_ref"] = NumericVector(hom.begin(), hom.end());
            ret["het"] = NumericVector(het.begin(), het.end());
            ret["hom_alt"] = NumericVector(alt.begin(), alt.end());
            return ret;
        }

        bool pass(int i, double min_maf, double max_maf, double min_cr, int min_mac, int max_mac, int m,
                  bool report_singletons) const {
            const double EPS = 1e-6;
            int left = 2 * hom[i] + het[i];
            int right = 2 * alt[i] + het[i];
            int sum = hom[i] + het[i] + alt[i];
            if (left < right) {
                std::swap(left, right);
            }
            double maf = (double)right / (left + right);
            double cr = (double)sum / m;
            bool singletons = report_singletons || (hom[i] + het[i] != 1 && het[i] + alt[i] != 1);
            return cr > min_cr && maf + EPS > min_maf && maf - EPS < max_maf && right >= min_mac &&
                      right <= max_mac && left > 0 && right > 0 && singletons;
        }

        int size() {
            return hom.size();
        }

        int get_hom(int i) const {
            return hom[i];
        }

        int get_het(int i) const {
            return het[i];
        }

        int get_alt(int i) const {
            return alt[i];
        }
    };

    int ceiling_devision(size_t x, size_t y) {
        return ((int)x + y - 1) / y;
    }

    class CountsReader {
        vector<size_t> samples;
        MemoryMappedScanner scanner;
        unsigned DP;
        unsigned GQ;

        int total_samples;
    public:
        CountsReader(const vector<size_t>& samples, const MemoryMappedScanner& scanner, unsigned DP, unsigned GQ, int tot_samples)
            :samples(samples), scanner(scanner), DP(DP), GQ(GQ), total_samples(tot_samples) {}

        CountsReader(const CountsReader&) = default;

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
            return std::tuple<int, int, int>{homref, het, hom};
        }
    };

    Counts parallel_read(const vector<size_t> &positions, const CountsReader& reader) {
        const int threads = 16;

        cxxpool::thread_pool pool{threads};
        std::vector<std::future<Counts>> futures;
        unsigned jobs_per_thread = (unsigned)ceiling_devision(positions.size(), threads);

        vector<size_t> thread_jobs;
        for (size_t i = 0; i < positions.size(); i++) {
            thread_jobs.push_back(positions[i]);
            if (thread_jobs.size() == jobs_per_thread || i == positions.size() - 1) {
                futures.push_back(pool.push([thread_jobs, reader]() -> Counts {
                    Counts counts = {};
                    for (size_t pos: thread_jobs) {
                        auto cts = reader.read(pos);
                        counts.push(std::get<0>(cts), std::get<1>(cts), std::get<2>(cts));
                    }
                    return counts;
                }));
                thread_jobs = vector<size_t>();
            }
        }

        Counts ret = {};

        for (auto& future: futures) {
            Counts done = future.get();
            ret.append(done);
        }
        return ret;
    }
}

namespace {
    List binaryScanResults(const vector<Variant>& variants, const unordered_map<Variant, bool>& requested,
            const vector<vcf::Range>& ranges, const Counts& counts, double min_maf, double max_maf, double cr,
            int min_mac, int max_mac, int n, bool report_singletons) {
        Counts cumulative;
        vector<string> names;
        vector<int> n_variants;
        int curr_range = 0;
        int range_entry = -1;
        for (size_t i = 0; i < variants.size(); i++) {
            if (i % 100 == 0) {
                Rcpp::checkUserInterrupt();
            }
            const Variant& var = variants[i];

            if (!counts.pass(i, min_maf, max_maf, cr, min_mac, max_mac, n, report_singletons)) {
                continue;
            }

            bool in_range = curr_range < (long)ranges.size();
            while (in_range && ranges[curr_range] < var.position()) {
                ++curr_range;
                range_entry = -1;
            }

            if (requested.find(var) != requested.end()) {
                bool as_is = requested.at(var);
                if (as_is) {
                    cumulative.push(counts.get_hom(i), counts.get_het(i), counts.get_alt(i));
                } else {
                    cumulative.push(counts.get_alt(i), counts.get_het(i), counts.get_hom(i));
                }
                names.push_back(as_is ? (string)var : (string)var.reversed());
                n_variants.push_back(1);
            }

            if (in_range && ranges[curr_range].includes(var.position())) {
                if (range_entry == -1) {
                    cumulative.push(0, 0, 0);
                    range_entry = cumulative.size() - 1;
                    names.push_back((string)ranges[curr_range]);
                    n_variants.push_back(0);
                }
                ++n_variants[range_entry];
                cumulative.add(range_entry, counts.get_hom(i), counts.get_het(i), counts.get_alt(i));
            }
        }
        List ret = cumulative.get_table();
        ret["names"] = CharacterVector(names.begin(), names.end());
        ret["n_variants"] = IntegerVector(n_variants.begin(), n_variants.end());
        return ret;
    }
}

// [[Rcpp::export]]
List parse_binary_file(const CharacterVector& variants, const CharacterVector& samples, const CharacterVector& regions,
        const CharacterVector& binary_file, const CharacterVector& metafile,
        const NumericVector& r_min_maf, const NumericVector& r_max_maf, const NumericVector& r_min_cr,
        const IntegerVector& r_min_mac, const IntegerVector& r_max_mac,
        const LogicalVector& report_singletons,
        const IntegerVector& requiredDP, const IntegerVector requiredGQ) {
    try {
        int DP = requiredDP[0];
        int GQ = requiredGQ[0];
        double min_maf = r_min_maf[0];
        double max_maf = r_max_maf[0];
        int min_mac = r_min_mac[0];
        int max_mac = r_max_mac[0];
        bool singletons = report_singletons[0];
        double cr = r_min_cr[0];

        RangeSet rangeSet;
        vector<vcf::Range> ranges;
        for (const char* s: regions) {
            auto r = vcf::Range::parseRange(std::string(s));
            ranges.push_back(r);
            rangeSet.insert(r);
        }

        MemoryMappedScanner scanner((string) binary_file[0]);
        ifstream fin(metafile[0]);

        string line;
        getline(fin, line);
        istringstream iss(line);
        unordered_map<string, int> sample_positions;
        size_t n = 0;
        while(iss) {
            string sample;
            if (!getline(iss, sample, '\t')) break;
            sample_positions.insert({sample, n++});
        }

        vector<size_t> positions;
        for (const char *s: samples) {
            if (sample_positions.find(s) == sample_positions.end()) {
                Rcpp::stop("Sample " + std::string(s) + " not found.");
            }
            positions.push_back(sample_positions[string(s)]);
        }

        unordered_map<Variant, bool> requested_variants;
        for (const char* var: variants) {
            auto variant = Variant::parseVariants(string(var))[0];
            requested_variants[variant] = true;
            requested_variants[variant.reversed()] = false;
        }

        vector<size_t> variant_pos;
        vector<Variant> found_variants;
        for (int i = 0; getline(fin, line); i++) {
            auto var = Variant::parseVariants(line)[0];
            if (rangeSet.includes(var.position()) || requested_variants.find(var) != requested_variants.end()) {
                found_variants.push_back(var);
                variant_pos.push_back(i);
            }
        }

        CountsReader reader(positions, scanner, DP, GQ, n);
        Counts counts = parallel_read(variant_pos, reader);

        int m = positions.size();
        List ret = binaryScanResults(found_variants, requested_variants, ranges, counts, min_maf, max_maf, cr,
                                     min_mac, max_mac, m, singletons);
        ret["total"] = positions.size();
        return ret;
    } catch (ParserException& e) {
        Rcpp::stop(e.get_message());
    }
}
