#include "vcf_parser.h"

#include <sstream>

namespace {
    using std::istream;
    using std::vector;
    using std::string;
}

namespace vcf {

    void VCFParser::registerHandler(std::unique_ptr<VariantsHandler>&& handler) {
        handlers.push_back(std::move(handler));
    }

    VCFParser::VCFParser(std::istream& input, const VCFFilter& filter) :input(input), filter(filter) {}

    std::vector<std::string> VCFParser::sample_names() {
        return samples;
    }

    void VCFParser::parseHeader() {
        int line_num = 0;
        string line;
        while(getline(input, line)) {
            ++line_num;
            if (line.substr(0, 2) == "##") {
                continue;
            }
            if (line.substr(0, 1) == "#") {
                std::istringstream in(line.substr(1));
                string token;
                for (int i = 0; getline(in, token, DELIM); i++) {
                    if (i < FIELDS.size()) {
                        if (token != FIELDS[i]) {
                            throw ParserException("Wrong header line: expected column " + FIELDS[i] +
                                                  "Found: " + token, line_num);
                        }
                    } else {
                        samples.push_back(token);
                    }
                }
            }
        }
    }
}
