#include "include/vcf_parser.h"

#include <algorithm>
#include <sstream>
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
        int first_allele, second_allele;
        std::istringstream iss(gt);
        iss >> first_allele;
        if (iss.eof()) {
            if (first_allele == 0) {
                return HOMREF;
            }
            return first_allele == allele ? HOM : MISSING;
        }
        char ch;
        iss >> ch;
        if (ch != DELIM_1 && ch != DELIM_2) {
            throw ParserException("Wrong GT format: " + gt);
        }
        iss >> second_allele;
        if (iss.fail()) {
            throw ParserException("Wrong GT format: " + gt);
        }
        return type(first_allele, second_allele, allele);
    }

    Allele Format::parse(const string& genotype, int allele, const VCFFilter& filter, VCFFilterStats& stats) {
        vector<string> parts = split(genotype, ':');
        try {
            string gt = parts[genotype_pos];
            if (gt == "." || gt == "./." || gt == ".|.") {
                stats.add(Stat::GT_MISS, 1);
                return {MISSING, 0, 0};
            }

            if (depth_pos >= parts.size() && qual_pos >= parts.size()) {
                throw ParserException("ignored");
            }

            int dp = depth_pos == -1 || parts[depth_pos] == "." ? 0 : stoi(parts[depth_pos]);
            int gq = qual_pos == -1 || parts[qual_pos] == "." ? 0 : stoi(parts[qual_pos]);

            if (!filter.apply(dp, gq)) {
                stats.add(Stat::DP_GQ, 1);
                return {MISSING, (unsigned)dp, (unsigned)gq};
            }
            Allele ret{parse_gt(gt, allele), (unsigned)dp, (unsigned)gq};
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
        auto tokens = split(*line, VCFParser::DELIM);

        if (tokens.size() != expected_ncols) {
            stats.add(Stat::WARNING, 1);
            throw ParserException("The row has " + std::to_string(tokens.size()) +
                                  " number of columns whereas header has " + std::to_string(expected_ncols));
        }
        Format format{tokens[FORMAT]};

        try {
            for (int sample : *indices) {
                alleles.push_back(format.parse(tokens.at(sample), variant + 1, *filter, stats));
            }
        } catch (ParserException& e) {
            corrupted = true;
            corruption_cause = e.get_message();
            throw e;
        }
    }

    std::vector<AlleleType> AlleleVector::vector() {
        resolve();
        std::vector<AlleleType> ret;
        ret.reserve(alleles.size());
        for (auto a: alleles) {
            ret.push_back(a.alleleType());
        }
        return ret;
    }

    void VCFParser::register_handler(std::shared_ptr<VariantsHandler> handler, int order) {
        handlers.emplace_back(handler, order);
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
        std::sort(handlers.begin(), handlers.end(), [](decltype(*handlers.begin())& l, decltype(*handlers.begin())& r){
            return l.second < r.second;
        });
        string line;
        std::shared_ptr<std::vector<size_t>> sample_indices = std::make_shared<std::vector<size_t>>(filtered_samples);
        std::shared_ptr<VCFFilter> vcf_filter = std::make_shared<VCFFilter>(filter);

        while (getline(input, line)) {
            Rcpp::checkUserInterrupt();
            ++line_num;
            if (line.empty() || std::all_of(line.begin(),line.end(),isspace)) {
                continue;
            }
            vector<string> tokens = split(line, DELIM, FIELDS.size());
            try {
                if (tokens.size() < FIELDS.size()) {
                    throw ParserException("The row is too short");
                }
                Position position = parse_position(tokens);
                vector<Variant> variants = parse_variants(tokens, position);
                stats.add(Stat::OVERALL, variants.size());

                if (tokens[FILTER] != "PASS") {
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
    }
}
