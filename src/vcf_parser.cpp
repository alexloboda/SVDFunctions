#include "vcf_parser.h"

#include <algorithm>
#include <sstream>
#include <iostream>

namespace {
    using namespace vcf;

    using std::istream;
    using std::vector;
    using std::string;
    using std::find;
    using std::pair;
    using std::stoi;

    vector<string> split(const string& line, char delim){
        vector<string> result;
        string curr;
        for (char ch : line) {
            if (ch != delim) {
                curr += ch;
            } else {
                if (curr.length() != 0) {
                    result.push_back(curr);
                    curr = "";
                }
            }
        }

        if (curr.length() != 0) {
            result.push_back(curr);
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

    class Format {
        const string DP_FIELD = "DP";
        const string GQ_FIELD = "GQ";
        const string GT_FIELD = "GT" ;
        const string AD_FIELD = "AD";

        const char DELIM_1 = '|';
        const char DELIM_2 = '/';

        const char MISSING_GT = '.';

        long depth_pos;
        long qual_pos;
        long genotype_pos;
        long ad_pos;

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

    public:
        Format(const string& format) {
            vector<string> parts = split(format, ':');
            find_pos(parts, DP_FIELD, depth_pos);
            find_pos(parts, GQ_FIELD, qual_pos);
            find_pos(parts, AD_FIELD, ad_pos);
            find_pos(parts, GT_FIELD, genotype_pos);
            if (genotype_pos == -1) {
                throw ParserException("No GT field available for a variant");
            }
        }

        AlleleType parse_gt(const string gt, int allele){
            if (find(gt.begin(), gt.end(), MISSING_GT) != gt.end()) {
                return MISSING;
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
                throw ParserException("Wrong GT format");
            }
            iss >> second_allele;
            if (iss.fail()) {
                throw ParserException("Wrong GT format");
            }
            return type(first_allele, second_allele, allele);
        }

        Allele parse(const string& genotype, int allele, const VCFFilter& filter) {
            vector<string> parts = split(genotype, ':');
            try {
                string gt = parts[genotype_pos];
                if (gt == "." || gt == "./.") {
                    return {MISSING, 0, 0};
                }
                if (ad_pos != -1) {
                    std::istringstream adstream(parts[ad_pos]);
                    int ref, alt;
                    char ch;
                    adstream >> ref >> ch >> alt;
                    if (!adstream.fail()) {
                        double ratio = ref / (double)(ref + alt);
                        if (ratio < 0.3 || ratio > 0.7) {
                            return {MISSING, 0, 0};
                        }
                    }
                }
                int dp = depth_pos == -1 ? 0 : stoi(parts[depth_pos]);
                int gq = qual_pos == -1 ? 0 : stoi(parts[qual_pos]);
                if (!filter.apply(dp, gq)) {
                    return {MISSING, (unsigned)dp, (unsigned)gq};
                }
                return {parse_gt(gt, allele), (unsigned)dp, (unsigned)gq};
            } catch (...) {
                throw ParserException("Wrong GT format");
            }
        }
    };
}

namespace vcf {

    void VCFParser::register_handler(std::shared_ptr<VariantsHandler> handler) {
        handlers.push_back(handler);
    }

    VCFParser::VCFParser(std::istream& input, const VCFFilter& filter) :input(input), filter(filter), line_num(0) {}

    std::vector<std::string> VCFParser::sample_names() {
        return samples;
    }

    vector<Variant> VCFParser::parse_variants(const vector<string>& tokens, const Position& position) {
        vector<Variant> variants;
        string ref = tokens[REF];
        vector<string> alts = split(tokens[ALT], ',');
        for (const string& alt: alts) {
            variants.emplace_back(position, ref, alt);
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
                vector<string> tokens = split(line, DELIM);
                for (int i = 0; i < tokens.size(); i++) {
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
    }

    void VCFParser::parse_genotypes() {
        string line;
        while (getline(input, line)) {
            ++line_num;
            vector<string> tokens = split(line, DELIM);
            if (tokens[FILTER] != "PASS") {
                continue;
            }
            try {
                Position position = parse_position(tokens);
                if (!filter.apply(position)) {
                    continue;
                }
                vector<Variant> variants = parse_variants(tokens, position);
                Format format(tokens[FORMAT]);
                for (int i = 0; i < variants.size(); i++) {
                    Variant& variant = variants[i];
                    if (filter.apply(variant)) {
                        vector<Allele> alleles;
                        for (int sample : filtered_samples) {
                            alleles.push_back(format.parse(tokens[sample], i + 1, filter));
                        }
                        for (auto& handler: handlers) {
                            handler->processVariant(variant, alleles);
                        }
                    }
                }
            } catch (const ParserException& e) {
                ParserException exception(e.get_message(), line_num);
                handle_error(exception);
            }
        }
    }

    void VCFParser::handle_error(const ParserException& e) {
        std::cerr << e.get_message() << std::endl;
    }
}
