#ifndef SRC_VCF_PARSER_H
#define SRC_VCF_PARSER_H

#include "vcf_primitives.h"
#include "vcf_filter.h"
#include "vcf_handlers.h"

namespace vcf {
    enum Field {
        CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
    };

    class VCFParser {
        static const char DELIM = '\t';
        const std::vector<std::string> FIELDS = {"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"};

        const VCFFilter& filter;

        std::vector<std::unique_ptr<VariantsHandler>> handlers;
        std::istream& input;
        std::vector<std::string> samples;

        int line_num;

    public:
        VCFParser(std::istream& input, const VCFFilter& filter);
        void parseHeader();
        void parseGenotypes();
        void registerHandler(std::unique_ptr<VariantsHandler>&& handler);

        std::vector<std::string> sample_names();
};

}


#endif //SRC_VCF_PARSER_H
