#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <future>
#include <vector>
#include <string>
#include <unordered_map>

#include "include/vcf_binary.h"
#include "include/matching.h"
#include "include/third-party/cxxpool.h"

using namespace Rcpp;
using namespace vcf;
using namespace std;

namespace {
    struct ScanFilters {
        double min_maf;
        double max_maf;
        double min_cr;
        int min_mac;
        int max_mac;
        bool report_singletons;
    };

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

        bool pass(int i, const ScanFilters& filters, int m) const {
            const double EPS = 1e-6;
            int left = 2 * hom[i] + het[i];
            int right = 2 * alt[i] + het[i];
            int sum = hom[i] + het[i] + alt[i];
            if (left < right) {
                std::swap(left, right);
            }
            double maf = (double)right / (left + right);
            double cr = (double)sum / m;
            bool singletons = filters.report_singletons || (hom[i] + het[i] != 1 && het[i] + alt[i] != 1);
            return cr > filters.min_cr && maf + EPS > filters.min_maf && maf - EPS < filters.max_maf &&
                      right >= filters.min_mac && right <= filters.max_mac && left > 0 && right > 0 && singletons;
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
    bool is_snv(const Variant& variant) {
        return variant.reference().size() == 1 && variant.alternative().size() == 1;
    }

    unordered_map<Variant, bool> build_requested_variants(const CharacterVector& variants) {
        unordered_map<Variant, bool> requested;
        for (const char* var: variants) {
            vector<Variant> parsed = Variant::parseVariants(string(var));
            if (parsed.empty()) {
                continue;
            }
            const Variant& variant = parsed[0];
            requested[variant] = true;
            if (is_snv(variant)) {
                requested[variant.reversed()] = false;
            }
        }
        return requested;
    }

    // Parses the requested scan regions into the ordered vector used to label
    // region rows and into the RangeSet used for fast membership tests.
    vector<vcf::Range> build_ranges(const CharacterVector& regions, RangeSet& rangeSet) {
        vector<vcf::Range> ranges;
        for (const char* s: regions) {
            auto r = vcf::Range::parseRange(std::string(s));
            ranges.push_back(r);
            rangeSet.insert(r);
        }
        return ranges;
    }

    // Reads the scalar quality-control thresholds passed from R into a single
    // ScanFilters bundle shared by both scanners.
    ScanFilters read_filters(const NumericVector& min_maf, const NumericVector& max_maf,
            const NumericVector& min_cr, const IntegerVector& min_mac, const IntegerVector& max_mac,
            const LogicalVector& report_singletons) {
        ScanFilters filters;
        filters.min_maf = min_maf[0];
        filters.max_maf = max_maf[0];
        filters.min_cr = min_cr[0];
        filters.min_mac = min_mac[0];
        filters.max_mac = max_mac[0];
        filters.report_singletons = report_singletons[0];
        return filters;
    }

    // Reads the tab-separated header line written by BinaryFileHandler and
    // returns the binary storage position of every sample.
    std::unordered_map<std::string, std::size_t> read_sample_positions(std::istream& meta, std::size_t& n_samples) {
        std::unordered_map<std::string, std::size_t> sample_positions;
        std::string header;
        std::getline(meta, header);
        std::istringstream iss(header);
        std::size_t n = 0;
        while (iss) {
            std::string sample;
            if (!std::getline(iss, sample, '\t')) {
                break;
            }
            sample_positions.emplace(sample, n++);
        }
        n_samples = n;
        return sample_positions;
    }

    // Walks the remaining variant lines of a metadata stream and records every
    // variant that is either explicitly requested or covered by a scan region,
    // together with its 0-based row index in the binary payload. Blank lines do
    // not map to a payload row and are skipped without advancing the index.
    void collect_scan_targets(std::istream& meta, RangeSet& rangeSet,
            const unordered_map<Variant, bool>& requested, vector<Variant>& found_variants,
            vector<size_t>& variant_pos) {
        std::string line;
        std::size_t row = 0;
        while (std::getline(meta, line)) {
            if (line.empty()) {
                continue;
            }
            auto var = Variant::parseVariants(line)[0];
            if (rangeSet.includes(var.position()) || requested.find(var) != requested.end()) {
                found_variants.push_back(var);
                variant_pos.push_back(row);
            }
            ++row;
        }
    }

    List binaryScanResults(const vector<Variant>& variants, const unordered_map<Variant, bool>& requested,
            const vector<vcf::Range>& ranges, const Counts& counts, const ScanFilters& filters, int n) {
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

            if (!counts.pass(i, filters, n)) {
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
        ScanFilters filters = read_filters(r_min_maf, r_max_maf, r_min_cr, r_min_mac, r_max_mac,
                                           report_singletons);

        RangeSet rangeSet;
        vector<vcf::Range> ranges = build_ranges(regions, rangeSet);

        MemoryMappedScanner scanner((string) binary_file[0]);
        ifstream fin(metafile[0]);

        size_t n = 0;
        auto sample_positions = read_sample_positions(fin, n);

        vector<size_t> positions;
        for (const char *s: samples) {
            auto it = sample_positions.find(string(s));
            if (it == sample_positions.end()) {
                Rcpp::stop("Sample " + std::string(s) + " not found.");
            }
            positions.push_back(it->second);
        }

        unordered_map<Variant, bool> requested_variants = build_requested_variants(variants);

        vector<size_t> variant_pos;
        vector<Variant> found_variants;
        collect_scan_targets(fin, rangeSet, requested_variants, found_variants, variant_pos);

        CountsReader reader(positions, scanner, DP, GQ, n);
        Counts counts = parallel_read(variant_pos, reader);

        int m = positions.size();
        List ret = binaryScanResults(found_variants, requested_variants, ranges, counts, filters, m);
        ret["total"] = positions.size();
        return ret;
    } catch (ParserException& e) {
        Rcpp::stop(e.get_message());
    }
}

namespace {
    constexpr std::size_t cluster_counts_bytes = sizeof(matching::ClusterCounts);
    constexpr int max_cluster_size = 255;
}

// [[Rcpp::export]]
List convert_to_cluster_binary(const CharacterVector& binary_file, const CharacterVector& metafile,
        const CharacterVector& samples, const IntegerVector& sample_clusters,
        const CharacterVector& cluster_labels, const CharacterVector& out_binary,
        const CharacterVector& out_meta, const IntegerVector& requiredDP, const IntegerVector& requiredGQ) {
    try {
        unsigned DP = requiredDP[0] < 0 ? 0u : (unsigned)requiredDP[0];
        unsigned GQ = requiredGQ[0] < 0 ? 0u : (unsigned)requiredGQ[0];
        std::size_t n_clusters = cluster_labels.size();
        if (n_clusters == 0) {
            Rcpp::stop("At least one cluster is required");
        }
        if (samples.size() != sample_clusters.size()) {
            Rcpp::stop("samples and sample_clusters must have the same length");
        }

        std::ifstream meta(metafile[0]);
        if (!meta) {
            Rcpp::stop("Cannot open metadata file");
        }
        std::size_t total_samples = 0;
        auto sample_positions = read_sample_positions(meta, total_samples);

        std::vector<int> position_cluster(total_samples, -1);
        std::vector<int> cluster_sizes(n_clusters, 0);
        for (int i = 0; i < samples.size(); i++) {
            std::string sample = (std::string)samples[i];
            auto it = sample_positions.find(sample);
            if (it == sample_positions.end()) {
                Rcpp::stop("Sample " + sample + " not found in metadata.");
            }
            int cluster = sample_clusters[i];
            if (cluster < 0 || (std::size_t)cluster >= n_clusters) {
                Rcpp::stop("Cluster index out of range");
            }
            position_cluster[it->second] = cluster;
            ++cluster_sizes[cluster];
        }
        for (std::size_t c = 0; c < n_clusters; c++) {
            if (cluster_sizes[c] > max_cluster_size) {
                Rcpp::stop("Each cluster must contain at most 255 samples");
            }
        }

        std::ifstream bin(binary_file[0], std::ios::binary);
        if (!bin) {
            Rcpp::stop("Cannot open binary file");
        }
        std::ofstream out_bin(out_binary[0], std::ios::binary | std::ios::trunc);
        if (!out_bin) {
            Rcpp::stop("Cannot open output binary file");
        }
        std::ofstream out_meta_stream(out_meta[0], std::ios::trunc);
        if (!out_meta_stream) {
            Rcpp::stop("Cannot open output metadata file");
        }

        for (std::size_t c = 0; c < n_clusters; c++) {
            out_meta_stream << (std::string)cluster_labels[c] << "\t";
        }
        out_meta_stream << "\n";
        for (std::size_t c = 0; c < n_clusters; c++) {
            out_meta_stream << cluster_sizes[c] << "\t";
        }
        out_meta_stream << "\n";

        std::vector<vcf::BinaryAllele> row(total_samples);
        std::vector<matching::ClusterCounts> block(n_clusters);
        const std::streamsize row_bytes = (std::streamsize)(total_samples * sizeof(vcf::BinaryAllele));
        std::string variant_line;
        std::size_t n_variants = 0;
        while (std::getline(meta, variant_line)) {
            if (variant_line.empty()) {
                continue;
            }
            bin.read(reinterpret_cast<char*>(row.data()), row_bytes);
            if (bin.gcount() != row_bytes) {
                Rcpp::stop("Binary file is truncated or inconsistent with metadata");
            }
            std::fill(block.begin(), block.end(), matching::ClusterCounts());
            if (n_variants % 100 == 0) {
                Rcpp::checkUserInterrupt();
            }
            for (std::size_t p = 0; p < total_samples; p++) {
                int cluster = position_cluster[p];
                if (cluster < 0) {
                    continue;
                }
                vcf::Allele allele = vcf::BinaryAllele::toAllele(row[p]);
                if (allele.DP() >= DP && allele.GQ() >= GQ) {
                    switch (allele.alleleType()) {
                        case vcf::HOMREF:
                            block[cluster][0] += 1;
                            break;
                        case vcf::HET:
                            block[cluster][1] += 1;
                            break;
                        case vcf::HOM:
                            block[cluster][2] += 1;
                            break;
                        default:
                            break;
                    }
                }
            }
            out_bin.write(reinterpret_cast<const char*>(block.data()),
                          (std::streamsize)(n_clusters * cluster_counts_bytes));
            out_meta_stream << variant_line << "\n";
            ++n_variants;
        }

        List ret;
        ret["variants"] = (int)n_variants;
        ret["clusters"] = (int)n_clusters;
        ret["samples"] = (int)samples.size();
        ret["total_samples"] = (int)total_samples;
        return ret;
    } catch (ParserException& e) {
        Rcpp::stop(e.get_message());
    }
}

// [[Rcpp::export]]
List parse_cluster_binary_file(const CharacterVector& variants, const CharacterVector& clusters,
        const CharacterVector& regions, const CharacterVector& binary_file, const CharacterVector& metafile,
        const NumericVector& r_min_maf, const NumericVector& r_max_maf, const NumericVector& r_min_cr,
        const IntegerVector& r_min_mac, const IntegerVector& r_max_mac, const LogicalVector& report_singletons) {
    try {
        ScanFilters filters = read_filters(r_min_maf, r_max_maf, r_min_cr, r_min_mac, r_max_mac,
                                           report_singletons);

        RangeSet rangeSet;
        vector<vcf::Range> ranges = build_ranges(regions, rangeSet);

        std::ifstream meta(metafile[0]);
        if (!meta) {
            Rcpp::stop("Cannot open metadata file");
        }

        std::string labels_line;
        std::getline(meta, labels_line);
        std::istringstream labels_iss(labels_line);
        std::unordered_map<std::string, std::size_t> cluster_positions;
        std::size_t n_clusters = 0;
        while (labels_iss) {
            std::string label;
            if (!std::getline(labels_iss, label, '\t')) {
                break;
            }
            cluster_positions.emplace(label, n_clusters++);
        }

        std::string sizes_line;
        if (!std::getline(meta, sizes_line)) {
            Rcpp::stop("Cluster metadata is missing the cluster sizes line");
        }
        std::istringstream sizes_iss(sizes_line);
        std::vector<int> cluster_sizes(n_clusters, 0);
        for (std::size_t c = 0; c < n_clusters; c++) {
            std::string size_token;
            if (!std::getline(sizes_iss, size_token, '\t') || size_token.empty()) {
                Rcpp::stop("Cluster metadata contains incomplete cluster sizes");
            }
            cluster_sizes[c] = std::stoi(size_token);
        }

        std::vector<std::size_t> selected;
        if (clusters.size() == 0) {
            selected.resize(n_clusters);
            std::iota(selected.begin(), selected.end(), 0);
        } else {
            for (const char* c: clusters) {
                auto it = cluster_positions.find(std::string(c));
                if (it == cluster_positions.end()) {
                    Rcpp::stop("Cluster " + std::string(c) + " not found.");
                }
                selected.push_back(it->second);
            }
        }

        long total = 0;
        for (std::size_t c: selected) {
            total += cluster_sizes[c];
        }

        unordered_map<Variant, bool> requested_variants = build_requested_variants(variants);

        vector<std::size_t> variant_pos;
        vector<Variant> found_variants;
        collect_scan_targets(meta, rangeSet, requested_variants, found_variants, variant_pos);

        std::ifstream bin(binary_file[0], std::ios::binary);
        if (!bin) {
            Rcpp::stop("Cannot open binary file");
        }
        const std::streamoff block_bytes = (std::streamoff)(n_clusters * cluster_counts_bytes);
        std::vector<matching::ClusterCounts> block(n_clusters);
        Counts counts;
        for (std::size_t k = 0; k < variant_pos.size(); k++) {
            if (k % 100 == 0) {
                Rcpp::checkUserInterrupt();
            }
            bin.seekg((std::streamoff)variant_pos[k] * block_bytes, std::ios::beg);
            bin.read(reinterpret_cast<char*>(block.data()), block_bytes);
            if (bin.gcount() != block_bytes) {
                Rcpp::stop("Cluster binary file is truncated or inconsistent with metadata");
            }
            int hom = 0, het = 0, alt = 0;
            for (std::size_t c: selected) {
                hom += block[c][0];
                het += block[c][1];
                alt += block[c][2];
            }
            counts.push(hom, het, alt);
        }

        List ret = binaryScanResults(found_variants, requested_variants, ranges, counts, filters, (int)total);
        ret["total"] = (int)total;
        return ret;
    } catch (ParserException& e) {
        Rcpp::stop(e.get_message());
    }
}
