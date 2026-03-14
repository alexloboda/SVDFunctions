#include "include/vcf_parser.h"
#include <Rcpp.h>
#include <boost/algorithm/string/predicate.hpp>
#include <cstdlib>
#include <fstream>
#include <chrono>
#include <cmath>
#include <memory>
#include <sstream>
#include <sys/stat.h>

#include "include/vcf_binary.h"
#include "include/third-party/zstr/zstr.hpp"
#include "include/third-party/zstr/strict_fstream.hpp"
#include "include/vcf_bgzf_parser.h"
#include "include/vcf_checkpoint.h"
#include "include/vcf_stats.h"
#include "include/vcf_predicting_handler.h"

namespace {
    using namespace Rcpp;
    using namespace vcf;
    using namespace std;
    using boost::algorithm::ends_with;

    constexpr int BGZF_COMPRESSION_TYPE = 2;

    static void hash_combine(uint64_t& seed, uint64_t value) {
        seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
    }

    static uint64_t stable_string_hash(const std::string& value) {
        uint64_t seed = 1469598103934665603ULL;
        for (unsigned char ch : value) {
            seed ^= static_cast<uint64_t>(ch);
            seed *= 1099511628211ULL;
        }
        return seed;
    }

    static uint64_t file_signature_bits(const std::string& path) {
        struct stat st {};
        if (::stat(path.c_str(), &st) != 0) {
            throw std::runtime_error("Failed to stat VCF file for checkpoint signature: " + path);
        }
        uint64_t seed = 0;
        hash_combine(seed, static_cast<uint64_t>(st.st_size));
        hash_combine(seed, static_cast<uint64_t>(st.st_mtime));
        return seed;
    }

    static uint64_t checkpoint_signature(const std::string& path,
                                         const CharacterVector& samples,
                                         const CharacterVector& bad_positions,
                                         const CharacterVector& variants,
                                         int dp,
                                         int gq,
                                         double missing_rate_threshold,
                                         unsigned int random_seed,
                                         int window_size,
                                         int rf_ntrees) {
        uint64_t seed = file_signature_bits(path);
        hash_combine(seed, stable_string_hash(path));
        hash_combine(seed, static_cast<uint64_t>(dp));
        hash_combine(seed, static_cast<uint64_t>(gq));
        hash_combine(seed, static_cast<uint64_t>(std::llround(missing_rate_threshold * 1e9)));
        hash_combine(seed, static_cast<uint64_t>(random_seed));
        hash_combine(seed, static_cast<uint64_t>(window_size));
        hash_combine(seed, static_cast<uint64_t>(rf_ntrees));
        for (auto sample : samples) {
            hash_combine(seed, stable_string_hash(Rcpp::as<std::string>(sample)));
        }
        for (auto position : bad_positions) {
            hash_combine(seed, stable_string_hash(Rcpp::as<std::string>(position)));
        }
        for (auto variant : variants) {
            hash_combine(seed, stable_string_hash(Rcpp::as<std::string>(variant)));
        }
        return seed;
    }

    static std::string checkpoint_prefix_from_dir(const std::string& checkpoint_dir) {
        return checkpoint_dir + "/predict_missing_checkpoint";
    }

    static bool checkpoint_crash_requested(const char* envvar, std::size_t generation) {
        const char* value = std::getenv(envvar);
        if (value == nullptr || *value == '\0') {
            return false;
        }
        try {
            return static_cast<std::size_t>(std::stoull(value)) == generation;
        } catch (...) {
            return false;
        }
    }

    static DataFrame make_loo_dataframe(const std::vector<ImputationLooRow>& rows) {
        CharacterVector variant(rows.size());
        IntegerVector n_observed(rows.size());
        IntegerVector n_missing(rows.size());
        NumericVector oob_mae(rows.size());
        NumericVector oob_rmse(rows.size());
        NumericVector rounded_acc(rows.size());
        for (size_t i = 0; i < rows.size(); i++) {
            variant[i] = rows[i].variant;
            n_observed[i] = static_cast<int>(rows[i].n_observed);
            n_missing[i] = static_cast<int>(rows[i].n_missing);
            oob_mae[i] = rows[i].oob_mae;
            oob_rmse[i] = rows[i].oob_rmse;
            rounded_acc[i] = rows[i].rounded_acc;
        }
        return DataFrame::create(
            _["variant"] = variant,
            _["n_observed"] = n_observed,
            _["n_missing"] = n_missing,
            _["oob_mae"] = oob_mae,
            _["oob_rmse"] = oob_rmse,
            _["rounded_acc"] = rounded_acc,
            _["stringsAsFactors"] = false
        );
    }

    class ProgressBar {
        bool enabled = false;
        std::shared_ptr<strict_fstream::ifstream> file;
        std::streamoff total_bytes = 0;
        Rcpp::Function set_txt_progress;
        Rcpp::Function close_fn;
        Rcpp::RObject bar;
        std::chrono::steady_clock::time_point last_update;
        double last_value = -1.0;
        std::string last_label;
    public:
                ProgressBar(bool enable, std::streamoff total, const std::shared_ptr<strict_fstream::ifstream>& file_in)
                        : enabled(enable && total > 0 && file_in),
                            file(file_in),
                            total_bytes(total),
                            set_txt_progress(Rcpp::Environment::base_env()["invisible"]),
                            close_fn(Rcpp::Environment::base_env()["invisible"]),
                            bar(R_NilValue) {
            if (!enabled) {
                return;
            }
            Rcpp::Environment base = Rcpp::Environment::base_env();
            Rcpp::Environment utils = Rcpp::Environment::namespace_env("utils");

            Rcpp::RObject stderr_obj = base["stderr"];
            Rcpp::RObject close_obj = base["close"];
            Rcpp::RObject txtpb_obj = utils["txtProgressBar"];
            Rcpp::RObject setpb_obj = utils["setTxtProgressBar"];

            if (stderr_obj.isNULL() || close_obj.isNULL() || txtpb_obj.isNULL() || setpb_obj.isNULL()) {
                enabled = false;
                return;
            }

            Rcpp::Function stderr_fn(stderr_obj);
            Rcpp::Function txtpb(txtpb_obj);
            set_txt_progress = Rcpp::Function(setpb_obj);
            close_fn = Rcpp::Function(close_obj);
            bar = txtpb(Rcpp::_["min"] = 0,
                        Rcpp::_["max"] = 100,
                        Rcpp::_["initial"] = 0,
                        Rcpp::_["title"] = "scanVCF",
                        Rcpp::_["style"] = 3,
                        Rcpp::_["file"] = stderr_fn());
            last_update = std::chrono::steady_clock::now();
        }

        void tick(const std::string& label) {
            if (!enabled) {
                return;
            }
            auto now = std::chrono::steady_clock::now();
            if (now - last_update < std::chrono::milliseconds(250)) {
                return;
            }
            last_update = now;
            std::streampos pos = file->tellg();
            if (pos < 0) {
                return;
            }
            double value = 100.0 * static_cast<double>(pos) / static_cast<double>(total_bytes);
            if (value > 100.0) {
                value = 100.0;
            }
            if (value < 0.0) {
                value = 0.0;
            }
            if (std::abs(value - last_value) < 0.1 && label == last_label) {
                return;
            }
            last_value = value;
            last_label = label;
            if (label.empty()) {
                set_txt_progress(bar, value);
            } else {
                set_txt_progress(bar, value, Rcpp::Named("label") = label);
            }
        }

        void finish() {
            if (!enabled) {
                return;
            }
            set_txt_progress(bar, 100.0);
            close_fn(bar);
            enabled = false;
        }

        ~ProgressBar() {
            if (!enabled) {
                return;
            }
            try {
                set_txt_progress(bar, 100.0);
                close_fn(bar);
            } catch (...) {
            }
        }
    };

    class Parser: public VCFParser {
        void handle_error(const vcf::ParserException& e) override {
            Rcpp::warning(e.get_message());
        }
    public:
        using VCFParser::VCFParser;
    };

    class RGenotypeMatrixHandler: public GenotypeMatrixHandler {
        bool checkpoint_enabled = false;
        bool checkpoint_complete = false;
        int64_t checkpoint_resume_offset = -1;
        int64_t first_data_offset = -1;
        std::size_t checkpoint_interval = 0;
        std::size_t loo_row_offset = 0;
        std::unique_ptr<PredictMissingCheckpoint> checkpoint;
        PredictMissingCheckpointManifest checkpoint_manifest;
    public:
        using GenotypeMatrixHandler::GenotypeMatrixHandler;

        void configure_checkpoint(const std::string& checkpoint_dir,
                                 std::size_t interval,
                                 uint64_t signature,
                                 bool resume_checkpoint,
                                 int compression,
                                 int64_t parser_first_data_offset) {
            if (compression != BGZF_COMPRESSION_TYPE) {
                throw std::runtime_error("checkpointDir requires a BGZF-compressed .vcf.gz input file");
            }

            checkpoint_enabled = true;
            checkpoint_interval = interval;
            first_data_offset = parser_first_data_offset;
            checkpoint.reset(new PredictMissingCheckpoint(checkpoint_prefix_from_dir(checkpoint_dir)));

            if (!resume_checkpoint) {
                checkpoint->reset();
            }

            checkpoint_manifest = PredictMissingCheckpointManifest();
            checkpoint_manifest.signature = signature;
            checkpoint_manifest.sample_count = samples.size();
            checkpoint_manifest.next_segment_id = 0;

            PredictMissingCheckpointManifest loaded_manifest;
            std::size_t visible_versions = 0;
            bool recovered_from_previous = false;
            std::string diagnostics;
            if (checkpoint->load_manifest(loaded_manifest, &visible_versions, &recovered_from_previous, &diagnostics)) {
                if (loaded_manifest.signature != signature) {
                    throw std::runtime_error("Existing checkpoint is incompatible with the current predictMissing request; rerun with resumeCheckpoint = FALSE or use a different checkpointDir");
                }
                if (loaded_manifest.sample_count != samples.size()) {
                    throw std::runtime_error("Existing checkpoint sample count does not match the current VCF header");
                }
                checkpoint_manifest = loaded_manifest;
                checkpoint_complete = loaded_manifest.complete;
                row_offset = loaded_manifest.flushed_rows;
                next_row_index = row_offset;
                loo_row_offset = loaded_manifest.flushed_loo_rows;
                checkpoint_resume_offset = loaded_manifest.flushed_rows == 0
                    ? first_data_offset
                    : loaded_manifest.resume_offset;
                if (!diagnostics.empty() && recovered_from_previous) {
                    Rcpp::warning(diagnostics);
                }
            } else {
                loo_row_offset = 0;
                checkpoint_resume_offset = first_data_offset;
            }
        }

        bool checkpoint_is_complete() const {
            return checkpoint_enabled && checkpoint_complete;
        }

        bool checkpoint_is_enabled() const {
            return checkpoint_enabled;
        }

        int64_t resume_offset() const {
            return checkpoint_resume_offset;
        }

        bool should_flush_checkpoint(const PredictingHandler& predicting_handler) const {
            return checkpoint_enabled && !checkpoint_complete &&
                predicting_handler.finalized_rows() >= checkpoint_manifest.flushed_rows + checkpoint_interval;
        }

        void flush_checkpoint(PredictingHandler* predicting_handler, bool complete) {
            if (!checkpoint_enabled || checkpoint_complete) {
                return;
            }

            std::vector<ImputationLooRow> empty_loo;
            const std::vector<ImputationLooRow>* loo_rows = &empty_loo;
            std::size_t stable_row_end = logical_size();
            int64_t next_resume_offset = -1;
            if (!complete && predicting_handler != nullptr) {
                predicting_handler->synchronize_loo();
                stable_row_end = predicting_handler->resume_row_index();
                next_resume_offset = predicting_handler->resume_offset();
                loo_rows = &predicting_handler->imputation_loo();
            } else if (!complete) {
                next_resume_offset = first_data_offset;
            }

            if (!complete && stable_row_end < checkpoint_manifest.flushed_rows) {
                throw std::runtime_error("Checkpoint resume row index moved backwards");
            }

            std::size_t stable_loo_end = checkpoint_manifest.flushed_loo_rows;
            std::size_t local_loo_index = stable_loo_end - loo_row_offset;
            while (local_loo_index < loo_rows->size() &&
                   (complete || (*loo_rows)[local_loo_index].logical_index < stable_row_end)) {
                ++stable_loo_end;
                ++local_loo_index;
            }

            PredictMissingCheckpointManifest next_manifest = checkpoint_manifest;
            next_manifest.generation += 1;

            if (stable_row_end > checkpoint_manifest.flushed_rows || stable_loo_end > checkpoint_manifest.flushed_loo_rows) {
                const std::size_t local_loo_from = checkpoint_manifest.flushed_loo_rows - loo_row_offset;
                const std::size_t local_loo_to = stable_loo_end - loo_row_offset;
                auto segment = checkpoint->write_segment(
                    next_manifest.next_segment_id,
                    variants,
                    gmatrix,
                    missing,
                    *loo_rows,
                    checkpoint_manifest.flushed_rows,
                    stable_row_end,
                    row_offset,
                    local_loo_from,
                    local_loo_to
                );
                next_manifest.segments.push_back(segment);
                next_manifest.next_segment_id += 1;
                next_manifest.flushed_rows = stable_row_end;
                next_manifest.flushed_loo_rows = stable_loo_end;

                if (checkpoint_crash_requested("SVDF_CHECKPOINT_TEST_ABORT_AFTER_SEGMENT_WRITE", next_manifest.generation)) {
                    throw std::runtime_error("Simulated checkpoint crash after segment write");
                }
            }

            next_manifest.complete = complete;
            next_manifest.resume_row_index = complete ? next_manifest.flushed_rows : stable_row_end;
            next_manifest.resume_offset = complete ? -1 : (next_manifest.flushed_rows == 0 ? first_data_offset : next_resume_offset);
            checkpoint->save_manifest(next_manifest);

            if (checkpoint_crash_requested("SVDF_CHECKPOINT_TEST_ABORT_AFTER_SAVE", next_manifest.generation)) {
                throw std::runtime_error("Simulated checkpoint crash after manifest save");
            }

            checkpoint_manifest = next_manifest;
            checkpoint_complete = complete;
            checkpoint_resume_offset = checkpoint_manifest.resume_offset;

            if (predicting_handler != nullptr && checkpoint_manifest.flushed_rows > row_offset) {
                predicting_handler->discard_checkpointed_prefix(checkpoint_manifest.flushed_rows);
            }
            if (predicting_handler != nullptr && checkpoint_manifest.flushed_loo_rows > loo_row_offset) {
                predicting_handler->discard_checkpointed_loo_prefix(checkpoint_manifest.flushed_loo_rows - loo_row_offset);
            }
            loo_row_offset = checkpoint_manifest.flushed_loo_rows;

            if (complete && predicting_handler != nullptr && !predicting_handler->imputation_loo().empty()) {
                PredictMissingCheckpointManifest tail_manifest = checkpoint_manifest;
                tail_manifest.generation += 1;
                auto tail_segment = checkpoint->write_segment(
                    tail_manifest.next_segment_id,
                    variants,
                    gmatrix,
                    missing,
                    predicting_handler->imputation_loo(),
                    checkpoint_manifest.flushed_rows,
                    checkpoint_manifest.flushed_rows,
                    row_offset,
                    0,
                    predicting_handler->imputation_loo().size()
                );
                tail_manifest.segments.push_back(tail_segment);
                tail_manifest.next_segment_id += 1;
                tail_manifest.flushed_loo_rows += predicting_handler->imputation_loo().size();
                checkpoint->replace_current_manifest(tail_manifest);
                predicting_handler->discard_checkpointed_loo_prefix(predicting_handler->imputation_loo().size());
                checkpoint_manifest = tail_manifest;
                loo_row_offset = checkpoint_manifest.flushed_loo_rows;
            }
        }

        std::vector<ImputationLooRow> load_checkpoint_loo() const {
            std::vector<ImputationLooRow> rows;
            if (checkpoint_enabled) {
                checkpoint->load_loo(rows);
            }
            return rows;
        }

        List result() {
            std::vector<std::string> row_names;
            std::vector<float> prefix_genotypes;
            std::vector<uint8_t> prefix_predicted;
            if (checkpoint_enabled) {
                checkpoint->load_prefix(samples.size(), row_names, prefix_genotypes, prefix_predicted);
            }

            const std::size_t prefix_rows = row_names.size();
            NumericMatrix res(prefix_rows + gmatrix.size(), samples.size());
            LogicalMatrix predicted(prefix_rows + missing.size(), samples.size());

            for (size_t i = 0; i < prefix_rows; i++) {
                for (size_t j = 0; j < samples.size(); j++) {
                    float val = prefix_genotypes[i * samples.size() + j];
                    if (val == to_int(vcf::MISSING)) {
                        val = NA_REAL;
                    }
                    res[j * (prefix_rows + gmatrix.size()) + i] = val;
                    predicted[j * (prefix_rows + missing.size()) + i] = prefix_predicted[i * samples.size() + j] != 0;
                }
            }

            for (size_t i = 0; i < gmatrix.size(); i++) {
                for (size_t j = 0; j < samples.size(); j++) {
                    float val = gmatrix[i][j];
                    if (val == to_int(vcf::MISSING)) {
                        val = NA_REAL;
                    }
                    res[j * (prefix_rows + gmatrix.size()) + prefix_rows + i] = val;
                    predicted[j * (prefix_rows + missing.size()) + prefix_rows + i] = missing[i][j];
                }
            }
            for_each(variants.begin(), variants.end(), [&row_names](Variant& v){
                row_names.push_back((string)v);
            });
            rownames(res) = CharacterVector(row_names.begin(), row_names.end());
            rownames(predicted) = CharacterVector(row_names.begin(), row_names.end());
            List ret;
            ret["genotype"] = res;
            ret["predicted"] = predicted;
            return ret;
        }
    };

    class RCallRateHandler: public CallRateHandler {
    public:
        using CallRateHandler::CallRateHandler;

        NumericMatrix result() {
            vector<string> non_empty;
            for (size_t i = 0; i < ranges.size(); i++) {
                if (n_variants[i] > 0) {
                    non_empty.push_back((std::string)ranges[i]);
                }
            }
            NumericMatrix result(non_empty.size(), samples.size());
            int curr = 0;
            for (size_t i = 0; i < ranges.size(); i++) {
                if (n_variants[i] == 0) {
                    continue;
                }
                for (size_t j = 0; j < samples.size(); j++) {
                    result[j * non_empty.size() + curr] = (double)call_rate_matrix[i][j] / n_variants[i];
                }
                ++curr;
            }
            rownames(result) = CharacterVector(non_empty.begin(), non_empty.end());
            return result;
        }
    };
}

VCFFilter filter(const CharacterVector& samples, const CharacterVector& bad_positions,
        int DP, int GQ) {
    VCFFilter filter(DP, GQ);

    if (samples.length() > 0) {
        vector<string> ss;
        for_each(samples.begin(), samples.end(), [&ss](const char *s) { ss.emplace_back(s); });
        filter.add_samples(ss);
    }

    if (bad_positions.length() > 0) {
        vector<Position> bads;
        for_each(bad_positions.begin(), bad_positions.end(), [&bads](const char *s) {
            bads.push_back(Position::parse_position(string(s)));
        });
        filter.add_bad_variants(bads);
    }

    return filter;
}

vector<vcf::Range> parse_regions(const CharacterVector& regions){
    vector<vcf::Range> ranges;
    for_each(regions.begin(), regions.end(), [&ranges](const char* str){
        ranges.push_back(vcf::Range::parseRange(string(str)));
    });
    return ranges;
}

// [[Rcpp::export]]
List parse_vcf(const CharacterVector& filename, const CharacterVector& samples,
               const CharacterVector& bad_positions, const CharacterVector& variants,
               const IntegerVector& DP, const IntegerVector& GQ, const LogicalVector& gmatrix,
               const LogicalVector& predictMissing, const CharacterVector& regions,
               const CharacterVector& binary_prefix, const NumericVector& missingRateThreshold,
               Rcpp::Nullable<int> seed, const IntegerVector& window_size,
               const IntegerVector& rf_ntrees, const CharacterVector& checkpoint_dir,
               const IntegerVector& checkpoint_interval, const LogicalVector& resume_checkpoint) {
    List ret;
    unsigned int random_seed = 42;
    if (seed.isNotNull()) {
        random_seed = (unsigned int)Rcpp::as<int>(seed);
    }
    int ws = 100;
    if (window_size.length() > 0) {
        ws = window_size[0];
    }
    if (ws < 3) {
        Rcpp::stop("window_size must be >= 3");
    }
    int ntrees = 50;
    if (rf_ntrees.length() > 0) {
        ntrees = rf_ntrees[0];
    }
    if (ntrees < 1) {
        Rcpp::stop("rf_ntrees must be >= 1");
    }
    std::string checkpoint_dir_value;
    if (checkpoint_dir.length() > 0) {
        checkpoint_dir_value = Rcpp::as<std::string>(checkpoint_dir[0]);
    }
    const std::size_t checkpoint_interval_value = checkpoint_interval.length() > 0
        ? static_cast<std::size_t>(checkpoint_interval[0])
        : 1000U;
    const bool resume_checkpoint_value = resume_checkpoint.length() == 0 || resume_checkpoint[0];
    try {
        const char *name = filename[0];
        VCFFilterStats stats;
        shared_ptr<RGenotypeMatrixHandler> gmatrix_handler;
        shared_ptr<BinaryFileHandler> binary_handler;
        shared_ptr<RCallRateHandler> callrate_handler;
        shared_ptr<PredictingHandler> predicting_handler;

        vector<Variant> vs;
        for_each(variants.begin(), variants.end(), [&vs](const char *s) {
            vector<Variant> parsed = Variant::parseVariants(string(s));
            vs.insert(vs.end(), parsed.begin(), parsed.end());
        });

        if (predictMissing[0]) {
            std::unique_ptr<BGZF, decltype(&bgzf_close)> bgzf_input(bgzf_open(name, "r"), &bgzf_close);
            if (!bgzf_input) {
                Rcpp::stop("Failed to open VCF file with BGZF reader");
            }

            BGZFVCFParser parser(bgzf_input.get(), filter(samples, bad_positions, DP[0], GQ[0]), stats,
                                 [](const vcf::ParserException& e) { Rcpp::warning(e.get_message()); });
            parser.parse_header();
            auto ss = parser.sample_names();
            parser.set_interrupt_every(2000);

            if (gmatrix[0]) {
                gmatrix_handler.reset(new RGenotypeMatrixHandler(ss, vs, stats, missingRateThreshold[0]));
                if (!checkpoint_dir_value.empty()) {
                    gmatrix_handler->configure_checkpoint(
                        checkpoint_dir_value,
                        checkpoint_interval_value,
                        checkpoint_signature(name, samples, bad_positions, variants, DP[0], GQ[0],
                                             missingRateThreshold[0], random_seed, ws, ntrees),
                        resume_checkpoint_value,
                        parser.compression(),
                        parser.first_data_offset()
                    );
                }
                parser.register_handler(gmatrix_handler, 1);
                if (!gmatrix_handler->checkpoint_is_complete()) {
                    predicting_handler = make_shared<PredictingHandler>(ss, *gmatrix_handler, 250000, ws,
                                                                        static_cast<std::size_t>(ntrees), random_seed);
                    parser.register_handler(predicting_handler, 2);
                }
            }

            if (!gmatrix_handler || !gmatrix_handler->checkpoint_is_complete()) {
                parser.parse_genotypes(gmatrix_handler && gmatrix_handler->checkpoint_is_enabled()
                                           ? gmatrix_handler->resume_offset()
                                           : -1,
                                       [&]() {
                                           if (gmatrix_handler && predicting_handler &&
                                               gmatrix_handler->should_flush_checkpoint(*predicting_handler)) {
                                               gmatrix_handler->flush_checkpoint(predicting_handler.get(), false);
                                           }
                                       });
                if (predicting_handler) {
                    predicting_handler->cleanup();
                }
                if (gmatrix_handler && gmatrix_handler->checkpoint_is_enabled()) {
                    gmatrix_handler->flush_checkpoint(predicting_handler.get(), true);
                }
            }

            ret["samples"] = CharacterVector(ss.begin(), ss.end());
            if (gmatrix[0]) {
                List geno = gmatrix_handler->result();
                if (predictMissing[0]) {
                    std::vector<ImputationLooRow> loo_rows = gmatrix_handler->checkpoint_is_enabled()
                        ? gmatrix_handler->load_checkpoint_loo()
                        : predicting_handler->imputation_loo();
                    geno["loo"] = make_loo_dataframe(loo_rows);
                }
                ret["genotype"] = geno;
            }
        } else {
            auto file = std::make_shared<strict_fstream::ifstream>(name, std::ios::in | std::ios::binary);
            file->seekg(0, std::ios::end);
            std::streamoff total_bytes = file->tellg();
            file->seekg(0, std::ios::beg);
            file->clear();

            auto zbuf = std::make_shared<zstr::istreambuf>(file->rdbuf());
            std::unique_ptr<std::istream> in(new std::istream(zbuf.get()));
            in->exceptions(std::ios_base::badbit);

            Parser parser(*in, filter(samples, bad_positions, DP[0], GQ[0]), stats);
            parser.parse_header();
            auto ss = parser.sample_names();

            Rcpp::Environment base = Rcpp::Environment::base_env();
            Rcpp::Function getOption = base["getOption"];
            Rcpp::Function interactive = base["interactive"];
            bool show_progress = Rcpp::as<bool>(getOption("svdf.progress", false));
            bool is_interactive = Rcpp::as<bool>(interactive());

            ProgressBar progress(show_progress && is_interactive, total_bytes, file);
            parser.set_progress_callback([&progress](const std::string& label) { progress.tick(label); }, 2000);
            parser.set_interrupt_every(2000);

            if (gmatrix[0]) {
                gmatrix_handler.reset(new RGenotypeMatrixHandler(ss, vs, stats, missingRateThreshold[0]));
                parser.register_handler(gmatrix_handler, 1);
            }

            if (regions.length() > 0) {
                callrate_handler.reset(new RCallRateHandler(ss, parse_regions(regions)));
                parser.register_handler(callrate_handler, 1);
            }

            if (binary_prefix.length() > 0) {
                string prefix = string(binary_prefix[0]);
                binary_handler.reset(new BinaryFileHandler(ss, prefix + "_bin", prefix + "_meta"));
                parser.register_handler(binary_handler, 1);
            }

            if (gmatrix_handler != nullptr || binary_handler != nullptr || callrate_handler != nullptr) {
                parser.parse_genotypes();
            }
            progress.finish();
            ret["samples"] = CharacterVector(ss.begin(), ss.end());
            if (gmatrix[0]) {
                ret["genotype"] = gmatrix_handler->result();
            }
            if (regions.length() > 0) {
                ret["callrate"] = callrate_handler->result();
            }
        }
        List ret_stats;
        for (Stat stat: vcf::statsList()) {
            ret_stats[to_string(stat)] = stats.value(stat);
        }
        ret["stats"] = ret_stats;
    } catch (ParserException& e) {
        Rcpp::stop(e.get_message());
    }
    return ret;
}

namespace {
    class Counts {
        vector<int> hom;
        vector<int> het;
        vector<int> alt;
    public:
        void push(int homc, int hetc, int altc) {
            hom.push_back(homc);
            het.push_back(hetc);
            alt.push_back(altc);
        }

        void add(int entry, int homc, int hetc, int altc) {
            hom[entry] += homc;
            het[entry] += hetc;
            alt[entry] += altc;
        }

        void append(const Counts& counts) {
            hom.insert(hom.end(), counts.hom.begin(), counts.hom.end());
            het.insert(het.end(), counts.het.begin(), counts.het.end());
            alt.insert(alt.end(), counts.alt.begin(), counts.alt.end());
        }

        List get_table() {
            List ret;
            ret["hom_ref"] = NumericVector(hom.begin(), hom.end());
            ret["het"] = NumericVector(het.begin(), het.end());
            ret["hom_alt"] = NumericVector(alt.begin(), alt.end());
            return ret;
        }

        bool pass(int i, double min_maf, double max_maf, double min_cr, int min_mac, int max_mac, int m,
                  bool report_singletons) const {
            const double EPS = 1e-6;
            int left = 2 * hom[i] + het[i];
            int right = 2 * alt[i] + het[i];
            int sum = hom[i] + het[i] + alt[i];
            if (left < right) {
                std::swap(left, right);
            }
            double maf = (double)right / (left + right);
            double cr = (double)sum / m;
            bool singletons = report_singletons || (hom[i] + het[i] != 1 && het[i] + alt[i] != 1);
            return cr > min_cr && maf + EPS > min_maf && maf - EPS < max_maf && right >= min_mac &&
                      right <= max_mac && left > 0 && right > 0 && singletons;
        }

        int size() {
            return hom.size();
        }

        int get_hom(int i) const {
            return hom[i];
        }

        int get_het(int i) const {
            return het[i];
        }

        int get_alt(int i) const {
            return alt[i];
        }
    };

    int ceiling_devision(size_t x, size_t y) {
        return ((int)x + y - 1) / y;
    }

    class CountsReader {
        vector<size_t> samples;
        MemoryMappedScanner scanner;
        unsigned DP;
        unsigned GQ;

        int total_samples;
    public:
        CountsReader(const vector<size_t>& samples, const MemoryMappedScanner& scanner, unsigned DP, unsigned GQ, int tot_samples)
            :samples(samples), scanner(scanner), DP(DP), GQ(GQ), total_samples(tot_samples) {}

        CountsReader(const CountsReader&) = default;

        std::tuple<int, int, int> read(size_t position) const {
            int homref = 0, het = 0, hom = 0;
            for (size_t sample_position: samples) {
                Allele allele = BinaryAllele::toAllele(scanner.scan(position * total_samples + sample_position));
                if (allele.DP() >= DP && allele.GQ() >= GQ) {
                    switch (allele.alleleType()) {
                        case HOM:
                            ++hom;
                            break;
                        case HET:
                            ++het;
                            break;
                        case HOMREF:
                            ++homref;
                            break;
                        default:
                            break;
                    }
                }
            }
            return std::tuple<int, int, int>{homref, het, hom};
        }
    };

    Counts parallel_read(const vector<size_t> &positions, const CountsReader& reader) {
        const int threads = 16;

        cxxpool::thread_pool pool{threads};
        std::vector<std::future<Counts>> futures;
        unsigned jobs_per_thread = (unsigned)ceiling_devision(positions.size(), threads);

        vector<size_t> thread_jobs;
        for (size_t i = 0; i < positions.size(); i++) {
            thread_jobs.push_back(positions[i]);
            if (thread_jobs.size() == jobs_per_thread || i == positions.size() - 1) {
                futures.push_back(pool.push([thread_jobs, reader]() -> Counts {
                    Counts counts = {};
                    for (size_t pos: thread_jobs) {
                        auto cts = reader.read(pos);
                        counts.push(std::get<0>(cts), std::get<1>(cts), std::get<2>(cts));
                    }
                    return counts;
                }));
                thread_jobs = vector<size_t>();
            }
        }

        Counts ret = {};

        for (auto& future: futures) {
            Counts done = future.get();
            ret.append(done);
        }
        return ret;
    }
}

namespace {
    List binaryScanResults(const vector<Variant>& variants, const unordered_map<Variant, bool>& requested,
            const vector<vcf::Range>& ranges, const Counts& counts, double min_maf, double max_maf, double cr,
            int min_mac, int max_mac, int n, bool report_singletons) {
        Counts cumulative;
        vector<string> names;
        vector<int> n_variants;
        int curr_range = 0;
        int range_entry = -1;
        for (size_t i = 0; i < variants.size(); i++) {
            if (i % 100 == 0) {
                Rcpp::checkUserInterrupt();
            }
            const Variant& var = variants[i];

            if (!counts.pass(i, min_maf, max_maf, cr, min_mac, max_mac, n, report_singletons)) {
                continue;
            }

            bool in_range = curr_range < (long)ranges.size();
            while (in_range && ranges[curr_range] < var.position()) {
                ++curr_range;
                range_entry = -1;
            }

            if (requested.find(var) != requested.end()) {
                bool as_is = requested.at(var);
                if (as_is) {
                    cumulative.push(counts.get_hom(i), counts.get_het(i), counts.get_alt(i));
                } else {
                    cumulative.push(counts.get_alt(i), counts.get_het(i), counts.get_hom(i));
                }
                names.push_back(as_is ? (string)var : (string)var.reversed());
                n_variants.push_back(1);
            }

            if (in_range && ranges[curr_range].includes(var.position())) {
                if (range_entry == -1) {
                    cumulative.push(0, 0, 0);
                    range_entry = cumulative.size() - 1;
                    names.push_back((string)ranges[curr_range]);
                    n_variants.push_back(0);
                }
                ++n_variants[range_entry];
                cumulative.add(range_entry, counts.get_hom(i), counts.get_het(i), counts.get_alt(i));
            }
        }
        List ret = cumulative.get_table();
        ret["names"] = CharacterVector(names.begin(), names.end());
        ret["n_variants"] = IntegerVector(n_variants.begin(), n_variants.end());
        return ret;
    }
}

// [[Rcpp::export]]
List parse_binary_file(const CharacterVector& variants, const CharacterVector& samples, const CharacterVector& regions,
        const CharacterVector& binary_file, const CharacterVector& metafile,
        const NumericVector& r_min_maf, const NumericVector& r_max_maf, const NumericVector& r_min_cr,
        const IntegerVector& r_min_mac, const IntegerVector& r_max_mac,
        const LogicalVector& report_singletons,
        const IntegerVector& requiredDP, const IntegerVector requiredGQ) {
    try {
        int DP = requiredDP[0];
        int GQ = requiredGQ[0];
        double min_maf = r_min_maf[0];
        double max_maf = r_max_maf[0];
        int min_mac = r_min_mac[0];
        int max_mac = r_max_mac[0];
        bool singletons = report_singletons[0];
        double cr = r_min_cr[0];

        RangeSet rangeSet;
        vector<vcf::Range> ranges;
        for (const char* s: regions) {
            auto r = vcf::Range::parseRange(std::string(s));
            ranges.push_back(r);
            rangeSet.insert(r);
        }

        MemoryMappedScanner scanner((string) binary_file[0]);
        ifstream fin(metafile[0]);

        string line;
        getline(fin, line);
        istringstream iss(line);
        unordered_map<string, int> sample_positions;
        size_t n = 0;
        while(iss) {
            string sample;
            if (!getline(iss, sample, '\t')) break;
            sample_positions.insert({sample, n++});
        }

        vector<size_t> positions;
        for (const char *s: samples) {
            if (sample_positions.find(s) == sample_positions.end()) {
                Rcpp::stop("Sample " + std::string(s) + " not found.");
            }
            positions.push_back(sample_positions[string(s)]);
        }

        unordered_map<Variant, bool> requested_variants;
        for (const char* var: variants) {
            auto variant = Variant::parseVariants(string(var))[0];
            requested_variants[variant] = true;
            requested_variants[variant.reversed()] = false;
        }

        vector<size_t> variant_pos;
        vector<Variant> found_variants;
        for (int i = 0; getline(fin, line); i++) {
            auto var = Variant::parseVariants(line)[0];
            if (rangeSet.includes(var.position()) || requested_variants.find(var) != requested_variants.end()) {
                found_variants.push_back(var);
                variant_pos.push_back(i);
            }
        }

        CountsReader reader(positions, scanner, DP, GQ, n);
        Counts counts = parallel_read(variant_pos, reader);

        int m = positions.size();
        List ret = binaryScanResults(found_variants, requested_variants, ranges, counts, min_maf, max_maf, cr,
                                     min_mac, max_mac, m, singletons);
        ret["total"] = positions.size();
        return ret;
    } catch (ParserException& e) {
        Rcpp::stop(e.get_message());
    }
}
