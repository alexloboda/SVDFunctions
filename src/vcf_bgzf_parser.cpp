#include "include/vcf_bgzf_parser.h"

#include <Rcpp.h>

#include <algorithm>
#include <cctype>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "htslib/kstring.h"

namespace {
    using namespace vcf;

    constexpr int BGZF_COMPRESSION_TYPE = 2;

    using std::string;
    using std::vector;

    vector<string> split_bgzf(const string& line, char delim, size_t max_num_tokens = 0) {
        size_t tokens = 1;
        for (char ch : line) {
            if (ch == delim) {
                ++tokens;
                if (tokens == max_num_tokens) {
                    break;
                }
            }
        }

        vector<string> result;
        result.reserve(tokens);
        size_t last = 0;
        for (size_t i = 0; i < line.length(); ++i) {
            if (line[i] == delim) {
                result.push_back(line.substr(last, i - last));
                last = i + 1;
                if (result.size() == max_num_tokens) {
                    return result;
                }
            }
        }
        if (last != line.length()) {
            result.push_back(line.substr(last));
        }
        return result;
    }

    Position parse_bgzf_position(const vector<string>& tokens) {
        Chromosome chr(tokens[CHROM]);
        int pos;
        try {
            pos = std::stoi(tokens[POS]);
        } catch (...) {
            throw ParserException("Can't read variant position");
        }
        return {chr, pos};
    }
}

namespace vcf {
    BGZFVCFParser::BGZFVCFParser(BGZF* input,
                                 const VCFFilter& filter,
                                 VCFFilterStats& stats,
                                 std::function<void(const ParserException&)> error_handler)
        : filter(filter),
          input(input),
          line_num(0),
          number_of_samples(0),
          stats(stats),
          data_offset(-1),
          error_handler(std::move(error_handler)),
          progress_every(0),
          interrupt_every(1000) {}

    void BGZFVCFParser::register_handler(std::shared_ptr<VariantsHandler> handler, int order) {
        handlers.emplace_back(std::move(handler), order);
    }

    void BGZFVCFParser::set_progress_callback(std::function<void(const std::string&)> callback,
                                              std::size_t every_n_lines) {
        progress_callback = std::move(callback);
        progress_every = every_n_lines;
    }

    void BGZFVCFParser::set_interrupt_every(std::size_t every_n_lines) {
        interrupt_every = every_n_lines;
    }

    std::vector<std::string> BGZFVCFParser::sample_names() const {
        return samples;
    }

    int64_t BGZFVCFParser::first_data_offset() const {
        return data_offset;
    }

    int BGZFVCFParser::compression() const {
        return bgzf_compression(input);
    }

    std::vector<Variant> BGZFVCFParser::parse_variants(const std::vector<std::string>& tokens,
                                                       const Position& position) {
        vector<Variant> variants;
        string ref = tokens[REF];
        vector<string> alts = split_bgzf(tokens[ALT], ',');
        for (const string& alt : alts) {
            Variant variant(position, ref, alt);
            if (is_of_interest(variant)) {
                variants.emplace_back(position, ref, alt);
            }
        }
        return variants;
    }

    bool BGZFVCFParser::is_of_interest(const Variant& var) {
        bool interesting = false;
        for (const auto& handler : handlers) {
            interesting |= handler.first->isOfInterest(var);
        }
        return interesting;
    }

    void BGZFVCFParser::parse_header() {
        kstring_t line = KS_INITIALIZE;
        while (bgzf_getline(input, '\n', &line) >= 0) {
            ++line_num;
            string current(line.s, line.l);
            if (current.substr(0, 2) == "##") {
                continue;
            }
            if (current.substr(0, 1) == "#") {
                current = current.substr(1);
                auto tokens = split_bgzf(current, VCFParser::DELIM);
                number_of_samples = static_cast<long>(tokens.size() - FIELDS.size());
                for (size_t i = 0; i < tokens.size(); ++i) {
                    const string& token = tokens[i];
                    if (i < FIELDS.size()) {
                        if (token != FIELDS[i]) {
                            ks_free(&line);
                            throw ParserException("Wrong header line: expected column " + FIELDS[i] + "Found: " + token,
                                                  line_num);
                        }
                    } else if (filter.apply(token)) {
                        samples.push_back(token);
                        filtered_samples.push_back(i);
                    }
                }
                data_offset = bgzf_tell(input);
                ks_free(&line);
                return;
            }
        }
        ks_free(&line);
        throw ParserException("No VCF header found in given file");
    }

    void BGZFVCFParser::parse_genotypes(int64_t resume_offset,
                                        const std::function<void()>& checkpoint_hook) {
        std::sort(handlers.begin(), handlers.end(), [](decltype(*handlers.cbegin()) l, decltype(*handlers.cbegin()) r) {
            return l.second < r.second;
        });

        const bool validating_resume = resume_offset >= 0 && resume_offset != data_offset;
        bool first_resumed_record_pending = validating_resume;

        if (resume_offset >= 0) {
            if (resume_offset != data_offset && bgzf_compression(input) != BGZF_COMPRESSION_TYPE) {
                throw ParserException("Checkpoint resume requires BGZF-compressed input");
            }
            if (resume_offset != data_offset && resume_offset < data_offset) {
                throw ParserException("Checkpoint resume offset precedes the first data row");
            }
            if (resume_offset != data_offset && bgzf_seek(input, resume_offset, SEEK_SET) < 0) {
                throw ParserException("Failed to seek BGZF stream to checkpoint offset");
            }
        }

        kstring_t line = KS_INITIALIZE;
        auto sample_indices = std::make_shared<std::vector<size_t>>(filtered_samples);
        auto vcf_filter = std::make_shared<VCFFilter>(filter);

        string last_position;
        if (progress_callback && progress_every > 0) {
            progress_callback(last_position);
        }

        while (true) {
            int64_t line_offset = bgzf_tell(input);
            int rc = bgzf_getline(input, '\n', &line);
            if (rc < 0) {
                if (first_resumed_record_pending) {
                    ks_free(&line);
                    throw ParserException("Checkpoint resume offset points past the last data row");
                }
                break;
            }

            ++line_num;
            if (interrupt_every > 0 && line_num % interrupt_every == 0) {
                Rcpp::checkUserInterrupt();
            }
            if (progress_callback && progress_every > 0 && line_num % progress_every == 0) {
                progress_callback(last_position);
            }

            string current(line.s, line.l);
            if (checkpoint_hook) {
                checkpoint_hook();
            }

            if (current.empty() || std::all_of(current.begin(), current.end(), [](unsigned char ch) {
                    return std::isspace(ch) != 0;
                })) {
                if (first_resumed_record_pending) {
                    ks_free(&line);
                    throw ParserException("Checkpoint resume offset does not point to a variant record");
                }
                continue;
            }

            if (first_resumed_record_pending && current[0] == '#') {
                ks_free(&line);
                throw ParserException("Checkpoint resume offset points into the VCF header instead of a variant record");
            }

            vector<string> tokens = split_bgzf(current, VCFParser::DELIM, FIELDS.size());
            try {
                if (tokens.size() < FIELDS.size()) {
                    throw ParserException("The row is too short");
                }
                Position position = parse_bgzf_position(tokens);
                last_position = tokens[CHROM] + ":" + tokens[POS];
                if (first_resumed_record_pending) {
                    first_resumed_record_pending = false;
                }
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

                auto line_pointer = std::make_shared<string>(std::move(current));
                for (size_t i = 0; i < variants.size(); ++i) {
                    Variant& variant = variants[i];
                    auto alleles = std::make_shared<AlleleVector>(line_pointer, sample_indices, vcf_filter, stats, i,
                                                                  FIELDS.size() + number_of_samples);
                    for (auto& handler : handlers) {
                        handler.first->processVariant(variant, alleles, line_offset);
                    }
                }
            } catch (const ParserException& e) {
                ParserException exception(e.get_message(), line_num);
                error_handler(exception);
            }
        }

        if (progress_callback && progress_every > 0) {
            progress_callback(last_position);
        }
        ks_free(&line);
    }
}