#ifndef SRC_VCF_PARSER_H
#define SRC_VCF_PARSER_H

#include "vcf_primitives.h"
#include "vcf_filter.h"

namespace vcf {
    class VCFParser {
        std::vector<std::string> fields;
        std::vector<Variant> variants_found;
        std::vector<std::string> samples_found;
        std::vector<std::vector<int>> gmatrix;
        std::unordered_map<std::string, int> fields_positions;

    public:
        explicit VCFParser(std::istream stream, VCFFilter& filter);
};

}


#endif //SRC_VCF_PARSER_H
