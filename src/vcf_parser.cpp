#include "include/vcf_parser.h"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <type_traits>
#include <Rcpp.h>

namespace {
    using namespace vcf;

    using std::istream;
    using std::vector;
    using std::string;
    using std::find;
    using std::pair;
    using std::stoi;

    // Rcpp conflict
    using vcf::MISSING;

    vector<std::string> split(const string& line, char delim, size_t max_num_tokens = 0){
        size_t tokens = 1;
        for (char ch: line) {
            if (ch == delim) {
                ++tokens;
                if (tokens == max_num_tokens) {
                    break;
                }
            }
        }

        vector<string> result;
        result.reserve(tokens);
        unsigned long last = 0;
        for (size_t i = 0; i < line.length(); i++) {
            char ch = line[i];
            if (ch == delim) {
                result.push_back(line.substr(last, i - last));
                last = i + 1;
                if (result.size() == max_num_tokens) {
                    return result;
                }
            }
        }

        if (last != line.length()) {
            result.push_back(line.substr(last, line.length() - last));
        }

        return result;
    }

    Position parse_position(const vector<string>& tokens) {
        Chromosome chr(tokens[CHROM]);
        int pos;
        try {
            pos = stoi(tokens[POS]);
        } catch (...) {
            throw ParserException("Can't read variant position");
        }
        return {chr, pos};
    }

    void find_pos(const vector<string>& tokens, const string& field, long& pos) {
        auto position = find(tokens.begin(), tokens.end(), field);
        if (position == tokens.end()) {
            pos = -1;
        } else {
            pos = position - tokens.begin();
        }
    }

    AlleleType type(int first, int second, int allele) {
        if (first > second) {
            std::swap(first, second);
        }
        if (first == second) {
            if (first == 0) {
                return HOMREF;
            } else if (first == allele) {
                return HOM;
            }
        } else {
            if (first == 0 && second == allele) {
                return HET;
            }
        }
        return MISSING;
    }

    struct Slice {
        size_t start;
        size_t end;
        bool found;
    };

    Slice locate_field(const string& value, long field_pos) {
        if (field_pos < 0) {
            return {0, 0, false};
        }

        size_t current = 0;
        size_t start = 0;
        for (size_t i = 0; i <= value.size(); i++) {
            if (i == value.size() || value[i] == ':') {
                if (current == (size_t)field_pos) {
                    return {start, i, true};
                }
                start = i + 1;
                current++;
            }
        }

        return {0, 0, false};
    }

    int parse_int_slice(const string& value, size_t start, size_t end, const string& original) {
        if (start >= end || (end - start == 1 && value[start] == '.')) {
            return 0;
        }

        int parsed = 0;
        for (size_t i = start; i < end; i++) {
            unsigned char ch = (unsigned char)value[i];
            if (!std::isdigit(ch)) {
                throw ParserException("Wrong GT format: " + original);
            }
            parsed = parsed * 10 + (value[i] - '0');
        }
        return parsed;
    }

    AlleleType parse_gt_slice(const string& gt, size_t start, size_t end, int allele) {
        if (start >= end) {
            throw ParserException("Wrong GT format: " + gt.substr(start, end - start));
        }

        auto parse_allele = [&](size_t& pos) {
            if (pos >= end || !std::isdigit((unsigned char)gt[pos])) {
                throw ParserException("Wrong GT format: " + gt.substr(start, end - start));
            }
            int parsed = 0;
            while (pos < end && std::isdigit((unsigned char)gt[pos])) {
                parsed = parsed * 10 + (gt[pos] - '0');
                pos++;
            }
            return parsed;
        };

        size_t pos = start;
        int first_allele = parse_allele(pos);
        if (pos == end) {
            if (first_allele == 0) {
                return HOMREF;
            }
            return first_allele == allele ? HOM : MISSING;
        }

        char sep = gt[pos++];
        if (sep != '|' && sep != '/') {
            throw ParserException("Wrong GT format: " + gt.substr(start, end - start));
        }

        int second_allele = parse_allele(pos);
        if (pos != end) {
            throw ParserException("Wrong GT format: " + gt.substr(start, end - start));
        }
        return type(first_allele, second_allele, allele);
    }
}

namespace vcf {
    Format::Format(const string& format) {
        vector<string> parts = split(format, ':');
        find_pos(parts, DP_FIELD, depth_pos);
        find_pos(parts, GQ_FIELD, qual_pos);
        find_pos(parts, AD_FIELD, ad_pos);
        find_pos(parts, GT_FIELD, genotype_pos);
        if (genotype_pos == -1) {
            throw ParserException("No GT field available for a variant");
        }
    }

    AlleleType Format::parse_gt(const string& gt, int allele){
        if (allele == 0) {
            return HOM;
        }
        return parse_gt_slice(gt, 0, gt.size(), allele);
    }

    Allele Format::parse(const string& genotype, int allele, const VCFFilter& filter, VCFFilterStats& stats) {
        return parse(genotype, 0, genotype.size(), allele, filter, stats);
    }

    Allele Format::parse(const string& genotype, size_t start, size_t end,
                         int allele, const VCFFilter& filter, VCFFilterStats& stats) {
        string original = genotype.substr(start, end - start);
        try {
            Slice gt = locate_field(original, genotype_pos);
            if (!gt.found) {
                throw ParserException("Wrong GT format: " + original);
            }

            if ((gt.end - gt.start == 1 && original[gt.start] == '.') ||
                (gt.end - gt.start == 3 && original[gt.start] == '.' && original[gt.start + 2] == '.' &&
                 (original[gt.start + 1] == '/' || original[gt.start + 1] == '|'))) {
                stats.add(Stat::GT_MISS, 1);
                return {MISSING, 0, 0};
            }

            Slice dp = locate_field(original, depth_pos);
            Slice gq = locate_field(original, qual_pos);

            if ((depth_pos >= 0 && !dp.found) || (qual_pos >= 0 && !gq.found)) {
                throw ParserException("ignored");
            }

            int dp_value = depth_pos == -1 ? 0 : parse_int_slice(original, dp.start, dp.end, original);
            int gq_value = qual_pos == -1 ? 0 : parse_int_slice(original, gq.start, gq.end, original);

            if (!filter.apply(dp_value, gq_value)) {
                stats.add(Stat::DP_GQ, 1);
                return {MISSING, (unsigned)dp_value, (unsigned)gq_value};
            }
            Allele ret{parse_gt_slice(original, gt.start, gt.end, allele), (unsigned)dp_value, (unsigned)gq_value};
            if (ret.alleleType() == HET) {
               // if (ad_pos != -1 && dp != 0) {
               //     std::istringstream adstream(parts[ad_pos]);
               //     int ref, alt;
               //     char ch;
               //     adstream >> ref;
               //     for (int i = 0; i < allele; i++) {
               //         while(!std::isdigit(adstream.peek())) {
               //             adstream >> ch;
               //         }
               //         adstream  >> alt;
               //     }
               //     if (!adstream.fail()) {
               //         double ref_ratio = ref / (double)dp;
               //         double alt_ratio = alt / (double)dp;
               //         if (ref_ratio < 0.3 || ref_ratio > 0.7) {
               //             stats.add(Stat::ALLELE_BALANCE, 1);
               //             return {MISSING, 0, 0};
               //         }
               //         if (alt_ratio < 0.3 || alt_ratio > 0.7) {
               //             stats.add(Stat::ALLELE_BALANCE, 1);
               //             return {MISSING, 0, 0};
               //         }
               //     }
               // }
            }
            return ret;
        } catch (...) {
            throw ParserException("Wrong GT format: " + genotype);
        }
    }

    AlleleVector::AlleleVector(std::shared_ptr<std::string>& line, std::shared_ptr<std::vector<size_t>>& indices,
            std::shared_ptr<vcf::VCFFilter>& filter, vcf::VCFFilterStats& stats, size_t variant, size_t ncols)
                                :line(line), indices(indices), filter(filter), stats(stats), variant(variant),
                                 expected_ncols(ncols) {}

    Allele AlleleVector::operator[](size_t i) {
        resolve();
        return alleles.at(i);
    }

    size_t AlleleVector::size() {
        resolve();
        return indices->size();
    }

    std::vector<Allele>::const_iterator AlleleVector::begin() {
        resolve();
        return alleles.begin();
    }

    std::vector<Allele>::const_iterator AlleleVector::end() {
        resolve();
        return alleles.end();
    }

    void AlleleVector::resolve() {
        if (corrupted) {
            throw ParserException(corruption_cause);
        }
        if (resolved) {
            return;
        }
        resolved = true;
        alleles.clear();
        alleles.reserve(indices->size());
        cached_unflipped.clear();
        cached_flipped.clear();
        cached_unflipped_ready = false;
        cached_flipped_ready = false;

        const string& raw = *line;
        std::vector<std::pair<size_t, size_t>> sample_bounds;
        sample_bounds.reserve(indices->size());
        size_t format_start = 0;
        size_t format_end = 0;
        bool format_found = false;
        size_t column = 0;
        size_t field_start = 0;
        size_t sample_cursor = 0;

        for (size_t i = 0; i <= raw.size(); i++) {
            if (i == raw.size() || raw[i] == VCFParser::DELIM) {
                if (column == FORMAT) {
                    format_start = field_start;
                    format_end = i;
                    format_found = true;
                }
                if (sample_cursor < indices->size() && (*indices)[sample_cursor] == column) {
                    sample_bounds.emplace_back(field_start, i);
                    sample_cursor++;
                }
                column++;
                field_start = i + 1;
            }
        }

        if (column != expected_ncols) {
            stats.add(Stat::WARNING, 1);
            throw ParserException("The row has " + std::to_string(column) +
                                  " number of columns whereas header has " + std::to_string(expected_ncols));
        }
        if (!format_found) {
            throw ParserException("No GT field available for a variant");
        }

        Format format{raw.substr(format_start, format_end - format_start)};

        try {
            for (const auto& bounds : sample_bounds) {
                alleles.push_back(format.parse(raw, bounds.first, bounds.second, variant + 1, *filter, stats));
            }
        } catch (ParserException& e) {
            corrupted = true;
            corruption_cause = e.get_message();
            throw e;
        }
    }

    std::vector<AlleleType> AlleleVector::vector(bool flipped) {
        resolve();

        if (!cached_unflipped_ready) {
            cached_unflipped.reserve(alleles.size());
            for (const auto& allele : alleles) {
                cached_unflipped.push_back(allele.alleleType());
            }
            cached_unflipped_ready = true;
        }

        if (!flipped) {
            return cached_unflipped;
        }

        if (!cached_flipped_ready) {
            cached_flipped = cached_unflipped;
            for (auto& allele : cached_flipped) {
                if (allele == HOM) {
                    allele = HOMREF;
                } else if (allele == HOMREF) {
                    allele = HOM;
                }
            }
            cached_flipped_ready = true;
        }

        return cached_flipped;
    }

    void VCFParser::register_handler(std::shared_ptr<VariantsHandler> handler, int order) {
        handlers.emplace_back(handler, order);
    }

    void VCFParser::set_progress_callback(std::function<void(const std::string&)> callback, std::size_t every_n_lines) {
        progress_callback = std::move(callback);
        progress_every = every_n_lines;
    }

    void VCFParser::set_interrupt_every(std::size_t every_n_lines) {
        interrupt_every = every_n_lines;
    }

    VCFParser::VCFParser(std::istream& input, const VCFFilter& filter, VCFFilterStats& stats) :filter(filter),
                         input(input), line_num(0), stats(stats){}

    std::vector<std::string> VCFParser::sample_names() {
        return samples;
    }

    vector<Variant> VCFParser::parse_variants(const vector<string>& tokens, const Position& position) {
        vector<Variant> variants;
        string ref = tokens[REF];
        vector<string> alts = split(tokens[ALT], ',');
        for (const string& alt: alts) {
            Variant variant(position, ref, alt);
            if (is_of_interest(variant)) {
                variants.emplace_back(position, ref, alt);
            }
        }
        return variants;
    }

    void VCFParser::parse_header() {
        string line;
        while (getline(input, line)) {
            ++line_num;
            if (line.substr(0, 2) == "##") {
                continue;
            }
            if (line.substr(0, 1) == "#") {
                line = line.substr(1);
                auto tokens = split(line, DELIM);
                number_of_samples = tokens.size() - FIELDS.size();
                for (size_t i = 0; i < tokens.size(); i++) {
                    const string& token = tokens[i];
                    if (i < FIELDS.size()) {
                        if (token != FIELDS[i]) {
                            throw ParserException("Wrong header line: expected column " + FIELDS[i] +
                                                  "Found: " + token, line_num);
                        }
                    } else {
                        if (filter.apply(token)) {
                            samples.push_back(token);
                            filtered_samples.push_back(i);
                        }
                    }
                }
                return;
            }
        }
        throw ParserException("No VCF header found in given file");
    }

    bool VCFParser::is_of_interest(const Variant& var) {
        bool interesting = false;
        for (const auto& handler: handlers) {
            interesting |= handler.first->isOfInterest(var);
        }
        return interesting;
    }

    void VCFParser::parse_genotypes() {
        std::sort(handlers.begin(), handlers.end(), [](decltype(*handlers.cbegin()) l, decltype(*handlers.cbegin()) r) {
            return l.second < r.second;
        });
        string line;
        std::shared_ptr<std::vector<size_t>> sample_indices = std::make_shared<std::vector<size_t>>(filtered_samples);
        std::shared_ptr<VCFFilter> vcf_filter = std::make_shared<VCFFilter>(filter);

        std::string last_position;
        if (progress_callback && progress_every > 0) {
            progress_callback(last_position);
        }
        while (getline(input, line)) {
            ++line_num;
            if (interrupt_every > 0 && line_num % interrupt_every == 0) {
                Rcpp::checkUserInterrupt();
            }
            if (progress_callback && progress_every > 0 && line_num % progress_every == 0) {
                progress_callback(last_position);
            }
            if (line.empty() || std::all_of(line.begin(),line.end(),isspace)) {
                continue;
            }
            vector<string> tokens = split(line, DELIM, FIELDS.size());
            try {
                if (tokens.size() < FIELDS.size()) {
                    throw ParserException("The row is too short");
                }
                Position position = parse_position(tokens);
                last_position = tokens[CHROM] + ":" + tokens[POS];
                vector<Variant> variants = parse_variants(tokens, position);
                stats.add(Stat::OVERALL, variants.size());

                if (tokens[FILTER] != "PASS" && tokens[FILTER] != ".") {
                    stats.add(Stat::NON_PASS, variants.size());
                    continue;
                }

                if (!filter.apply(position)) {
                    stats.add(Stat::BANNED, variants.size());
                    continue;
                }

                if (variants.empty()) {
                    continue;
                }
                std::shared_ptr<std::string> line_pointer = std::make_shared<std::string>(std::move(line));

                for (size_t i = 0; i < variants.size(); i++) {
                    Variant& variant = variants[i];
                    auto alleles = std::make_shared<AlleleVector>(line_pointer, sample_indices, vcf_filter, stats, i,
                                                                  FIELDS.size() + number_of_samples);
                    for (auto& handler: handlers) {
                        handler.first->processVariant(variant, alleles);
                    }
                }
            } catch (const ParserException& e) {
                ParserException exception(e.get_message(), line_num);
                handle_error(exception);
            }
        }

        if (progress_callback && progress_every > 0) {
            progress_callback(last_position);
        }
    }
}
