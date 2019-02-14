#include "vcf_parser.h"

#include <sstream>
#include <iostream>

namespace {
    using namespace vcf;

    using std::istream;
    using std::vector;
    using std::string;

    vector<string> split(const string& line, char delim){
        vector<string> result;
        string curr;
        for (int i = 0; i < line.length(); i++) {
            if (line[i] != delim) {
                curr += line[i];
            } else {
                if (curr.length() != 0) {
                    result.push_back(curr);
                    curr = "";
                }
            }
        }
        return result;
    }
}

namespace vcf {

    void VCFParser::registerHandler(std::unique_ptr<VariantsHandler>&& handler) {
        handlers.push_back(std::move(handler));
    }

    VCFParser::VCFParser(std::istream& input, const VCFFilter& filter) :input(input), filter(filter), line_num(0) {}

    std::vector<std::string> VCFParser::sample_names() {
        return samples;
    }

    vector<Variant> VCFParser::parse_variants(const vector<string>& tokens) {
        Chromosome chr(tokens[CHROM]);
        int pos;
        try {
            pos = std::stoi(tokens[POS]);
        } catch (...) {
            throw ParserException("Can't read variant position");
        }
        Position position(chr, pos);
        vector<Variant> variants;
        string ref = tokens[REF];
        vector<string> alts = split(tokens[ALT], ',');
        for (const string& alt: alts) {
            Variant variant(position, ref, alt);
            if (filter.apply(variant)) {
                variants.emplace_back(position, ref, alt);
            }
        }
        return variants;
    }

    void VCFParser::parseHeader() {
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
            }
        }
    }

    void VCFParser::parseGenotypes() {
        string line;
        while (getline(input, line)) {
            ++line_num;
            vector<string> tokens = split(line, DELIM);
            if (tokens[QUAL] != "PASS") {
                continue;
            }
            try {
                vector<Variant> variants = parse_variants(tokens);
            } catch (const ParserException& e) {
                ParserException exception(e.get_message(), line_num);
                std::cerr << e.get_message() << std::endl;
            }
        }
    }
}
