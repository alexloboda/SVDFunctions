#ifndef SRC_VCF_PARSER_H
#define SRC_VCF_PARSER_H

#include "vcf_primitives.h"
#include "vcf_filter.h"
#include "vcf_handlers.h"
#include "vcf_stats.h"

namespace vcf {
    enum Field {
        CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT
    };

    class VCFParser {
        static const char DELIM = '\t';
        const std::vector<std::string> FIELDS = {"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"};

        VCFFilter filter;

        std::vector<std::pair<std::shared_ptr<VariantsHandler>, int>> handlers;
        std::istream& input;
        std::vector<std::string> samples;
        std::vector<int> filtered_samples;

        int line_num;
        long number_of_samples;

        VCFFilterStats& stats;

        std::vector<Variant> parse_variants(const std::vector<std::string>& tokens, const Position& position);
        virtual void handle_error(const ParserException& e) = 0;
        bool is_of_interest(const Variant& var);

    public:
        VCFParser(std::istream& input, const VCFFilter& filter, VCFFilterStats& stats);
        void parse_header();
        void parse_genotypes();
        void register_handler(std::shared_ptr<VariantsHandler> handler, int order);

        std::vector<std::string> sample_names();
};

}


#endif //SRC_VCF_PARSER_H
