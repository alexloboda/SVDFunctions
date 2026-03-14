#ifndef SRC_VCF_BGZF_PARSER_H
#define SRC_VCF_BGZF_PARSER_H

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "htslib/bgzf.h"

#include "vcf_filter.h"
#include "vcf_handlers.h"
#include "vcf_parser.h"
#include "vcf_stats.h"

namespace vcf {
    class BGZFVCFParser {
        VCFFilter filter;
        BGZF* input;
        std::vector<std::pair<std::shared_ptr<VariantsHandler>, int>> handlers;
        std::vector<std::string> samples;
        std::vector<std::size_t> filtered_samples;
        int line_num;
        long number_of_samples;
        VCFFilterStats& stats;
        int64_t data_offset;
        std::function<void(const ParserException&)> error_handler;
        std::function<void(const std::string&)> progress_callback;
        std::size_t progress_every;
        std::size_t interrupt_every;

        std::vector<Variant> parse_variants(const std::vector<std::string>& tokens, const Position& position);
        bool is_of_interest(const Variant& var);
    public:
        BGZFVCFParser(BGZF* input,
                      const VCFFilter& filter,
                      VCFFilterStats& stats,
                      std::function<void(const ParserException&)> error_handler);

        void parse_header();
        void parse_genotypes(int64_t resume_offset = -1,
                             const std::function<void()>& checkpoint_hook = {});
        void register_handler(std::shared_ptr<VariantsHandler> handler, int order);
        void set_progress_callback(std::function<void(const std::string&)> callback,
                                   std::size_t every_n_lines = 10000);
        void set_interrupt_every(std::size_t every_n_lines);

        std::vector<std::string> sample_names() const;
        int64_t first_data_offset() const;
        int compression() const;
    };
}

#endif //SRC_VCF_BGZF_PARSER_H