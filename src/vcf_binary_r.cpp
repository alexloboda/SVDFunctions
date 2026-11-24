#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <atomic>
#include <charconv>
#include <chrono>
#include <iterator>
#include <numeric>
#include <optional>
#include <ranges>
#include <span>
#include <string_view>
#include <system_error>
#include <thread>
#include <tuple>
#include <type_traits>
#include <future>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>

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

            // The bound is rechecked on every step: a variant past the last
            // range moves curr_range to the end of the vector.
            while (curr_range < (long)ranges.size() && ranges[curr_range] < var.position()) {
                ++curr_range;
                range_entry = -1;
            }
            bool in_range = curr_range < (long)ranges.size();

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

// The per-cluster scanner. Both of its passes, over the metadata lines and over
// the payload rows, are cut into tasks by the size of the input alone and the
// task results are put together in file order. The number of threads therefore
// only decides how many tasks run at once, never what the scan returns.
//
// Both files are memory mapped. The mio sources own the mappings and release
// them when they go out of scope; everything else only looks into them through
// std::string_view and std::span, which own nothing.
namespace {
    constexpr std::size_t METADATA_CHUNK_BYTES = 256u << 10;
    constexpr std::size_t PAYLOAD_ROWS_PER_TASK = 4096;

    // Maps a whole file read-only; `what` names the file in the error message.
    // An empty file cannot be mapped, so it is left as an empty source instead
    // of being an error.
    template<typename Source>
    Source map_file(const std::string& path, const std::string& what) {
        std::ifstream probe(path, std::ios::binary | std::ios::ate);
        if (!probe) {
            Rcpp::stop("Cannot open " + what);
        }
        Source source;
        if (probe.tellg() > 0) {
            std::error_code error;
            source.map(path, error);
            if (error) {
                Rcpp::stop("Cannot map " + what + ": " + error.message());
            }
        }
        return source;
    }

    // Hands out the lines of a text one after another, the way std::getline does
    // for a stream, but as views into the text instead of copies. The metadata
    // is written in text mode, so a file made on Windows ends its lines with
    // "\r\n" and the carriage return goes along with the line feed.
    class LineReader {
        std::string_view rest;
    public:
        explicit LineReader(std::string_view text) : rest(text) {}

        // The next line, or nothing once the text is used up.
        std::optional<std::string_view> next() {
            if (rest.empty()) {
                return std::nullopt;
            }
            std::string_view line = rest.substr(0, rest.find('\n'));
            rest.remove_prefix(std::min(line.size() + 1, rest.size()));
            if (line.ends_with('\r')) {
                line.remove_suffix(1);
            }
            return line;
        }

        // The text that has not been handed out yet.
        std::string_view remaining() const { return rest; }
    };

    // Runs task(0) ... task(n_tasks - 1) on up to `threads` workers and returns
    // the results in task order. Tasks run off the main thread, so they must not
    // touch the R API; user interrupts are polled here instead. The pool drains
    // its queue before it joins, hence the flag that turns whatever is still
    // queued into no-ops when this function is left early.
    template<typename Task>
    auto run_tasks(std::size_t n_tasks, int threads, const Task& task) {
        using Result = std::invoke_result_t<const Task&, std::size_t>;
        std::vector<Result> results;
        results.reserve(n_tasks);

        std::size_t workers = threads > 0 ? static_cast<std::size_t>(threads)
                                          : std::thread::hardware_concurrency();
        workers = std::max<std::size_t>(1, std::min(workers, n_tasks));
        if (workers == 1) {
            for (std::size_t i = 0; i < n_tasks; i++) {
                Rcpp::checkUserInterrupt();
                results.push_back(task(i));
            }
            return results;
        }

        struct RaiseOnExit {
            std::atomic<bool>& flag;
            ~RaiseOnExit() { flag = true; }
        };
        std::atomic<bool> cancelled{false};
        cxxpool::thread_pool pool(workers);
        RaiseOnExit cancel{cancelled};

        std::vector<std::future<Result>> futures;
        futures.reserve(n_tasks);
        for (std::size_t i = 0; i < n_tasks; i++) {
            futures.push_back(pool.push([&task, &cancelled, i]() {
                return cancelled ? Result() : task(i);
            }));
        }
        for (auto& future: futures) {
            while (future.wait_for(std::chrono::milliseconds(100)) != std::future_status::ready) {
                Rcpp::checkUserInterrupt();
            }
            results.push_back(future.get());
        }
        return results;
    }

    // What operator>> takes for the gap between two fields of a line.
    bool is_space(char c) {
        return c == ' ' || (c >= '\t' && c <= '\r');
    }

    // The part of a text before its first whitespace.
    std::string_view first_field(std::string_view text) {
        return text.substr(0, std::ranges::find_if(text, is_space) - text.begin());
    }

    std::string_view skip_spaces(std::string_view text) {
        text.remove_prefix(std::ranges::find_if_not(text, is_space) - text.begin());
        return text;
    }

    // Drops a leading "chr", in any letter case, the way Chromosome does.
    std::string_view without_chr(std::string_view text) {
        constexpr std::string_view prefix = "chr";
        auto same_letter = [](char c, char lower) { return c == lower || c == lower - ('a' - 'A'); };
        if (std::ranges::equal(text.substr(0, prefix.size()), prefix, same_letter)) {
            text.remove_prefix(prefix.size());
        }
        return text;
    }

    // A number read off the front of a text, and the text that follows it.
    struct LeadingNumber {
        int value;
        std::string_view rest;
    };

    // Reads the decimal digits a text starts with. Nothing when it does not
    // start with a digit, or when the number does not fit an int.
    std::optional<LeadingNumber> leading_number(std::string_view text) {
        // std::from_chars would also take a minus sign.
        if (text.empty() || text.front() < '0' || text.front() > '9') {
            return std::nullopt;
        }
        int value = 0;
        auto [stop, error] = std::from_chars(text.data(), text.data() + text.size(), value);
        if (error != std::errc()) {
            return std::nullopt;
        }
        return LeadingNumber{value, text.substr(stop - text.data())};
    }

    // The position a variant line starts with, and the text that follows it.
    struct LeadingPosition {
        Position position;
        std::string_view rest;
    };

    // Reads the position off the front of a variant line when it is spelled the
    // plain way: "chr#:#" with an autosome number. Nothing for any other
    // spelling, which is then left to Variant::parseVariants.
    std::optional<LeadingPosition> leading_position(std::string_view line) {
        auto chromosome = leading_number(without_chr(line));
        if (!chromosome || chromosome->value < 1 || chromosome->value > 22 ||
                !chromosome->rest.starts_with(':')) {
            return std::nullopt;
        }
        auto coordinate = leading_number(chromosome->rest.substr(1));
        if (!coordinate || (!coordinate->rest.empty() && !is_space(coordinate->rest.front()))) {
            return std::nullopt;
        }
        return LeadingPosition{Position(Chromosome(chromosome->value), coordinate->value), coordinate->rest};
    }

    // Parses a whole variant line the way the sequential scanner does, so both
    // scanners accept and reject exactly the same lines.
    Variant parse_variant_line(std::string_view line) {
        std::string text(line);
        vector<Variant> parsed = Variant::parseVariants(text);
        if (parsed.empty()) {
            throw ParserException("Malformed variant in metadata: " + text);
        }
        return parsed.front();
    }

    // One variant line of the metadata, read in place.
    //
    // Most lines of a large file are not scan targets, and the position alone
    // tells. So a plainly spelled position is read on its own, and the rest of
    // the line is only looked at when the variant is asked for. A line that
    // spells its position any other way is parsed in full straight away.
    class VariantLine {
        std::string_view line;
        std::optional<LeadingPosition> plain;
        std::optional<Variant> parsed;
    public:
        // Throws ParserException when the line had to be parsed in full and is malformed.
        explicit VariantLine(std::string_view line) : line(line), plain(leading_position(line)) {
            if (!plain) {
                parsed = parse_variant_line(line);
            }
        }

        Position position() const {
            return plain ? plain->position : parsed->position();
        }

        // Throws ParserException when the alleles of the line are malformed.
        Variant variant() const {
            if (parsed) {
                return *parsed;
            }
            // REF, then the first of the comma separated ALTs, which is the
            // variant Variant::parseVariants puts first.
            std::string_view fields = skip_spaces(plain->rest);
            std::string_view ref = first_field(fields);
            std::string_view alts = skip_spaces(fields.substr(ref.size()));
            if (ref.empty() || alts.empty()) {
                return parse_variant_line(line);
            }
            std::string_view alt = alts.substr(0, alts.find(','));
            return {plain->position, std::string(ref), std::string(alt)};
        }
    };

    struct ScanTargets {
        vector<Variant> variants;
        // Payload rows of the variants, counted from the first line scanned.
        vector<std::size_t> rows;
        // Number of payload rows the scanned text spans, i.e. its non-blank lines.
        std::size_t n_rows = 0;
    };

    // Scans a chunk of metadata lines for the variants that are requested
    // explicitly or covered by a scan region.
    ScanTargets scan_metadata_chunk(std::string_view chunk, const RangeSet& ranges,
            const unordered_map<Variant, bool>& requested,
            const std::unordered_set<Position>& requested_positions) {
        ScanTargets targets;
        LineReader lines(chunk);
        while (std::optional<std::string_view> line = lines.next()) {
            if (line->empty()) {
                continue;
            }
            std::size_t row = targets.n_rows++;

            VariantLine entry(*line);
            Position position = entry.position();
            bool in_range = ranges.includes(position);
            if (!in_range && !requested_positions.contains(position)) {
                continue;
            }
            Variant variant = entry.variant();
            if (in_range || requested.contains(variant)) {
                targets.variants.push_back(std::move(variant));
                targets.rows.push_back(row);
            }
        }
        return targets;
    }

    // Parallel counterpart of the stream based collect_scan_targets: `text`
    // holds the variant lines of the mapped metadata, which are cut into chunks
    // on line boundaries. A chunk cannot know how many rows precede it, so it
    // numbers its rows from zero and the offsets are added while the chunks are
    // put together.
    ScanTargets collect_scan_targets(std::string_view text, const RangeSet& ranges,
            const unordered_map<Variant, bool>& requested, int threads) {
        std::unordered_set<Position> requested_positions;
        for (const Variant& variant: requested | std::views::keys) {
            requested_positions.insert(variant.position());
        }

        vector<std::string_view> chunks;
        while (text.size() > METADATA_CHUNK_BYTES) {
            std::size_t line_feed = text.find('\n', METADATA_CHUNK_BYTES);
            if (line_feed == std::string_view::npos) {
                break;
            }
            chunks.push_back(text.substr(0, line_feed + 1));
            text.remove_prefix(line_feed + 1);
        }
        if (!text.empty()) {
            chunks.push_back(text);
        }

        auto parts = run_tasks(chunks.size(), threads, [&](std::size_t i) {
            return scan_metadata_chunk(chunks[i], ranges, requested, requested_positions);
        });

        std::size_t total = 0;
        for (const ScanTargets& part: parts) {
            total += part.rows.size();
        }
        ScanTargets targets;
        targets.variants.reserve(total);
        targets.rows.reserve(total);
        for (ScanTargets& part: parts) {
            std::ranges::move(part.variants, std::back_inserter(targets.variants));
            for (std::size_t row: part.rows) {
                targets.rows.push_back(targets.n_rows + row);
            }
            targets.n_rows += part.n_rows;
            part = ScanTargets();
        }
        return targets;
    }

    // Sums the selected clusters of a payload row: the counterpart of
    // CountsReader for the per-cluster binary. A row holds one
    // hom-ref/het/hom-alt byte triple per cluster. The payload is mapped once
    // and read by all threads.
    class ClusterCountsReader {
        mio::ummap_source payload;
        // Where the triple of every selected cluster starts within a row.
        vector<std::size_t> offsets;
        std::size_t row_bytes;
    public:
        ClusterCountsReader(const std::string& path, const vector<std::size_t>& clusters,
                std::size_t n_clusters)
            : payload(map_file<mio::ummap_source>(path, "binary file")),
              row_bytes(n_clusters * cluster_counts_bytes) {
            offsets.reserve(clusters.size());
            for (std::size_t cluster: clusters) {
                offsets.push_back(cluster * cluster_counts_bytes);
            }
        }

        // Whether the file is long enough to hold the given row.
        bool holds(std::size_t row) const {
            return (row + 1) * row_bytes <= payload.size();
        }

        std::tuple<int, int, int> read(std::size_t row) const {
            std::span<const unsigned char> bytes(payload.data(), payload.size());
            std::span<const unsigned char> counts = bytes.subspan(row * row_bytes, row_bytes);
            int hom_ref = 0, het = 0, hom_alt = 0;
            for (std::size_t offset: offsets) {
                hom_ref += counts[offset];
                het += counts[offset + 1];
                hom_alt += counts[offset + 2];
            }
            return {hom_ref, het, hom_alt};
        }
    };

    // Reads the requested payload rows. They are ascending, as they were
    // collected in file order, so every task walks its own stretch of the file
    // front to back and the last row decides whether the file holds them all.
    Counts read_cluster_counts(const ClusterCountsReader& reader, const vector<std::size_t>& rows,
            int threads) {
        if (!rows.empty() && !reader.holds(rows.back())) {
            Rcpp::stop("Cluster binary file is truncated or inconsistent with metadata");
        }

        std::size_t n_tasks = (rows.size() + PAYLOAD_ROWS_PER_TASK - 1) / PAYLOAD_ROWS_PER_TASK;
        auto parts = run_tasks(n_tasks, threads, [&](std::size_t task) {
            std::size_t first = task * PAYLOAD_ROWS_PER_TASK;
            std::size_t count = std::min(PAYLOAD_ROWS_PER_TASK, rows.size() - first);
            Counts counts;
            for (std::size_t row: std::span(rows).subspan(first, count)) {
                auto [hom_ref, het, hom_alt] = reader.read(row);
                counts.push(hom_ref, het, hom_alt);
            }
            return counts;
        });

        Counts counts;
        for (const Counts& part: parts) {
            counts.append(part);
        }
        return counts;
    }
}

// [[Rcpp::export]]
List parse_cluster_binary_file(const CharacterVector& variants, const CharacterVector& clusters,
        const CharacterVector& regions, const CharacterVector& binary_file, const CharacterVector& metafile,
        const NumericVector& r_min_maf, const NumericVector& r_max_maf, const NumericVector& r_min_cr,
        const IntegerVector& r_min_mac, const IntegerVector& r_max_mac, const LogicalVector& report_singletons,
        int threads) {
    try {
        ScanFilters filters = read_filters(r_min_maf, r_max_maf, r_min_cr, r_min_mac, r_max_mac,
                                           report_singletons);

        RangeSet rangeSet;
        vector<vcf::Range> ranges = build_ranges(regions, rangeSet);

        mio::mmap_source meta = map_file<mio::mmap_source>((std::string)metafile[0], "metadata file");
        LineReader meta_lines{std::string_view(meta.data(), meta.size())};

        std::string labels_line(meta_lines.next().value_or(""));
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

        std::optional<std::string_view> sizes_line = meta_lines.next();
        if (!sizes_line) {
            Rcpp::stop("Cluster metadata is missing the cluster sizes line");
        }
        std::istringstream sizes_iss{std::string(*sizes_line)};
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

        ScanTargets targets = collect_scan_targets(meta_lines.remaining(), rangeSet, requested_variants,
                                                   threads);

        ClusterCountsReader reader((std::string)binary_file[0], selected, n_clusters);
        Counts counts = read_cluster_counts(reader, targets.rows, threads);

        List ret = binaryScanResults(targets.variants, requested_variants, ranges, counts, filters, (int)total);
        ret["total"] = (int)total;
        return ret;
    } catch (ParserException& e) {
        Rcpp::stop(e.get_message());
    }
}
