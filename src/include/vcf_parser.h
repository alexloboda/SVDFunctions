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

    const static std::vector<std::string> FIELDS = {"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"};

    class Format {
        const std::string DP_FIELD = "DP";
        const std::string GQ_FIELD = "GQ";
        const std::string GT_FIELD = "GT" ;
        const std::string AD_FIELD = "AD";

        const char DELIM_1 = '|';
        const char DELIM_2 = '/';

        long depth_pos;
        long qual_pos;
        long genotype_pos;
        long ad_pos;

    public:
        Format(const std::string& format);
        AlleleType parse_gt(const std::string& gt, int allele);
        Allele parse(const std::string& genotype, int allele, const VCFFilter& filter, VCFFilterStats& stats);
    };

    class AlleleVector {
        std::shared_ptr<std::string> line;
        std::shared_ptr<std::vector<size_t>> indices;
        std::shared_ptr<VCFFilter> filter;
        vcf::VCFFilterStats& stats;
        std::vector<Allele> alleles;
        size_t variant;
        bool resolved = false;

        void resolve();
    public:
        AlleleVector(std::shared_ptr<std::string>& line, std::shared_ptr<std::vector<size_t>>& indices,
                std::shared_ptr<VCFFilter>& filter, VCFFilterStats& stats, size_t variant);
        AlleleVector(const AlleleVector&) = delete;

        std::vector<Allele>::const_iterator begin();
        std::vector<Allele>::const_iterator end();
        size_t size();
        Allele operator[](size_t i);
    };

    class VCFParser {
        VCFFilter filter;

        std::vector<std::pair<std::shared_ptr<VariantsHandler>, int>> handlers;
        std::istream& input;
        std::vector<std::string> samples;
        std::vector<size_t> filtered_samples;

        int line_num;
        long number_of_samples;

        VCFFilterStats& stats;

        std::vector<Variant> parse_variants(const std::vector<std::string>& tokens, const Position& position);
        virtual void handle_error(const ParserException& e) = 0;
        bool is_of_interest(const Variant& var);

    public:
        static const char DELIM = '\t';

        VCFParser(std::istream& input, const VCFFilter& filter, VCFFilterStats& stats);
        void parse_header();
        void parse_genotypes();
        void register_handler(std::shared_ptr<VariantsHandler> handler, int order);

        std::vector<std::string> sample_names();
};

}


#endif //SRC_VCF_PARSER_H
