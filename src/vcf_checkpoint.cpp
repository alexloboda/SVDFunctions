#include "include/vcf_checkpoint.h"

#include <cerrno>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#else
#include <fcntl.h>
#include <unistd.h>
#endif

#include "include/vcf_predicting_handler.h"

namespace {
    using std::size_t;
    using std::string;

    string make_path(const string& prefix, const char* suffix) {
        return prefix + suffix;
    }

    bool local_file_exists(const string& path) {
        std::ifstream input(path, std::ios::binary);
        return input.good();
    }

    uint64_t local_file_size(const string& path) {
        std::ifstream input(path, std::ios::binary | std::ios::ate);
        if (!input.good()) {
            throw std::runtime_error("Checkpoint file is not readable: " + path);
        }
        return static_cast<uint64_t>(input.tellg());
    }

    size_t count_lines(const string& path) {
        std::ifstream input(path);
        if (!input.good()) {
            throw std::runtime_error("Checkpoint file is not readable: " + path);
        }
        size_t count = 0;
        string line;
        while (std::getline(input, line)) {
            ++count;
        }
        if (!input.eof()) {
            throw std::runtime_error("Failed while reading checkpoint file: " + path);
        }
        return count;
    }

    string make_segment_path(const string& prefix, size_t segment_id, const char* suffix) {
        return prefix + "_predict_segment_" + std::to_string(segment_id) + suffix;
    }

    vcf::PredictMissingCheckpointManifest parse_manifest_text(const string& text) {
        vcf::PredictMissingCheckpointManifest loaded;
        std::istringstream input(text);
        string line;
        bool magic_ok = false;
        while (std::getline(input, line)) {
            std::istringstream iss(line);
            string key;
            string value;
            if (!std::getline(iss, key, '\t') || !std::getline(iss, value)) {
                continue;
            }
            if (key == "magic") {
                if (value == "SVDFunctionsPredictMissingCheckpointV3") {
                    magic_ok = true;
                } else if (value == "SVDFunctionsPredictMissingCheckpointV2" ||
                           value == "SVDFunctionsPredictMissingCheckpointV1") {
                    throw std::runtime_error("Legacy checkpoint format is incompatible with the current predictMissing output; restart with resumeCheckpoint = FALSE");
                }
            } else if (key == "generation") {
                loaded.generation = static_cast<size_t>(std::stoull(value));
            } else if (key == "signature") {
                loaded.signature = static_cast<uint64_t>(std::stoull(value));
            } else if (key == "sample_count") {
                loaded.sample_count = static_cast<size_t>(std::stoull(value));
            } else if (key == "flushed_rows") {
                loaded.flushed_rows = static_cast<size_t>(std::stoull(value));
            } else if (key == "flushed_loo_rows") {
                loaded.flushed_loo_rows = static_cast<size_t>(std::stoull(value));
            } else if (key == "resume_offset") {
                loaded.resume_offset = static_cast<int64_t>(std::stoll(value));
            } else if (key == "resume_row_index") {
                loaded.resume_row_index = static_cast<size_t>(std::stoull(value));
            } else if (key == "complete") {
                loaded.complete = (value == "1");
            } else if (key == "next_segment_id") {
                loaded.next_segment_id = static_cast<size_t>(std::stoull(value));
            } else if (key == "segment") {
                vcf::PredictMissingCheckpointSegment segment;
                std::istringstream segment_stream(value);
                string rows;
                string loo_rows;
                if (!std::getline(segment_stream, value, '\t') ||
                    !std::getline(segment_stream, rows, '\t') ||
                    !std::getline(segment_stream, loo_rows)) {
                    throw std::runtime_error("Checkpoint manifest segment entry is malformed");
                }
                segment.id = static_cast<size_t>(std::stoull(value));
                segment.rows = static_cast<size_t>(std::stoull(rows));
                segment.loo_rows = static_cast<size_t>(std::stoull(loo_rows));
                loaded.segments.push_back(segment);
            }
        }

        if (!magic_ok) {
            throw std::runtime_error("Checkpoint manifest is invalid or has unsupported format");
        }
        if (loaded.generation == 0) {
            loaded.generation = 1;
        }
        return loaded;
    }

    string serialize_manifest_text(const vcf::PredictMissingCheckpointManifest& manifest) {
        std::ostringstream output;
        output << "magic\tSVDFunctionsPredictMissingCheckpointV3\n";
        output << "generation\t" << manifest.generation << "\n";
        output << "signature\t" << manifest.signature << "\n";
        output << "sample_count\t" << manifest.sample_count << "\n";
        output << "flushed_rows\t" << manifest.flushed_rows << "\n";
        output << "flushed_loo_rows\t" << manifest.flushed_loo_rows << "\n";
        output << "resume_offset\t" << manifest.resume_offset << "\n";
        output << "resume_row_index\t" << manifest.resume_row_index << "\n";
        output << "complete\t" << (manifest.complete ? 1 : 0) << "\n";
        output << "next_segment_id\t" << manifest.next_segment_id << "\n";
        for (const auto& segment : manifest.segments) {
            output << "segment\t" << segment.id << "\t" << segment.rows << "\t" << segment.loo_rows << "\n";
        }
        return output.str();
    }

    string parent_dir(const string& path) {
        const auto pos = path.find_last_of("/\\");
        if (pos == string::npos) {
            return ".";
        }
        return path.substr(0, pos);
    }

    [[noreturn]] void throw_io_error(const string& action, const string& path) {
        throw std::runtime_error(action + ": " + path + ": " + std::strerror(errno));
    }

    void sync_file(const string& path) {
#ifdef _WIN32
        int fd = _open(path.c_str(), _O_BINARY | _O_RDWR);
        if (fd < 0) {
            throw_io_error("Failed to open checkpoint file for sync", path);
        }
        if (_commit(fd) != 0) {
            const int saved_errno = errno;
            _close(fd);
            errno = saved_errno;
            throw_io_error("Failed to sync checkpoint file", path);
        }
        if (_close(fd) != 0) {
            throw_io_error("Failed to close checkpoint file after sync", path);
        }
#else
        int fd = ::open(path.c_str(), O_RDONLY);
        if (fd < 0) {
            throw_io_error("Failed to open checkpoint file for sync", path);
        }
        if (::fsync(fd) != 0) {
            const int saved_errno = errno;
            ::close(fd);
            errno = saved_errno;
            throw_io_error("Failed to sync checkpoint file", path);
        }
        if (::close(fd) != 0) {
            throw_io_error("Failed to close checkpoint file after sync", path);
        }
#endif
    }

    void sync_directory(const string& path) {
#ifndef _WIN32
        string dir = parent_dir(path);
        int flags = O_RDONLY;
#ifdef O_DIRECTORY
        flags |= O_DIRECTORY;
#endif
        int fd = ::open(dir.c_str(), flags);
        if (fd < 0) {
            throw_io_error("Failed to open checkpoint directory for sync", dir);
        }
        if (::fsync(fd) != 0) {
            const int saved_errno = errno;
            ::close(fd);
            errno = saved_errno;
            throw_io_error("Failed to sync checkpoint directory", dir);
        }
        if (::close(fd) != 0) {
            throw_io_error("Failed to close checkpoint directory after sync", dir);
        }
#else
        (void) path;
#endif
    }

    void remove_if_exists(const string& path) {
        std::remove(path.c_str());
    }

    void write_manifest_atomically(const string& manifest_tmp_file,
                                   const string& manifest_file,
                                   const string& manifest_previous_file,
                                   const vcf::PredictMissingCheckpointManifest& manifest,
                                   bool rotate_previous) {
        std::ofstream output(manifest_tmp_file, std::ios::trunc);
        if (!output.good()) {
            throw std::runtime_error("Failed to write checkpoint manifest: " + manifest_tmp_file);
        }
        output << serialize_manifest_text(manifest);
        output.flush();
        if (!output.good()) {
            throw std::runtime_error("Failed to flush checkpoint manifest: " + manifest_tmp_file);
        }
        output.close();
        if (!output.good()) {
            throw std::runtime_error("Failed to close checkpoint manifest: " + manifest_tmp_file);
        }
        sync_file(manifest_tmp_file);

        if (rotate_previous && local_file_exists(manifest_file)) {
            if (std::rename(manifest_file.c_str(), manifest_previous_file.c_str()) != 0) {
                throw std::runtime_error("Failed to rotate checkpoint manifest: " + string(std::strerror(errno)));
            }
        } else if (!rotate_previous && local_file_exists(manifest_file)) {
            remove_if_exists(manifest_file);
        }

        if (std::rename(manifest_tmp_file.c_str(), manifest_file.c_str()) != 0) {
            throw std::runtime_error("Failed to publish checkpoint manifest: " + string(std::strerror(errno)));
        }
        sync_directory(manifest_file);
    }

    string describe_manifest(const string& slot_name, const vcf::PredictMissingCheckpointManifest& manifest) {
        std::ostringstream out;
        out << slot_name << " generation " << manifest.generation
            << " (rows=" << manifest.flushed_rows
            << ", loo_rows=" << manifest.flushed_loo_rows
            << ", complete=" << (manifest.complete ? "yes" : "no") << ")";
        return out.str();
    }

    void validate_manifest(const vcf::PredictMissingCheckpointManifest& manifest, const string& prefix) {
        size_t total_rows = 0;
        size_t total_loo_rows = 0;
        size_t max_segment_id = 0;
        bool have_segments = false;
        for (const auto& segment : manifest.segments) {
            total_rows += segment.rows;
            total_loo_rows += segment.loo_rows;
            max_segment_id = have_segments ? std::max(max_segment_id, segment.id) : segment.id;
            have_segments = true;
        }

        if (manifest.sample_count == 0) {
            throw std::runtime_error("Checkpoint manifest is missing sample_count for prefix: " + prefix);
        }
        if (total_rows != manifest.flushed_rows) {
            throw std::runtime_error("Checkpoint manifest row count does not match committed segments for prefix: " + prefix);
        }
        if (total_loo_rows != manifest.flushed_loo_rows) {
            throw std::runtime_error("Checkpoint manifest loo row count does not match committed segments for prefix: " + prefix);
        }
        if (!manifest.complete && manifest.flushed_rows > 0 && manifest.resume_offset < 0) {
            throw std::runtime_error("Checkpoint manifest is missing a valid resume offset for prefix: " + prefix);
        }
        if (!manifest.complete && manifest.resume_row_index > manifest.flushed_rows) {
            throw std::runtime_error("Checkpoint manifest has an invalid resume row index for prefix: " + prefix);
        }
        if (manifest.complete && manifest.resume_offset != -1) {
            throw std::runtime_error("Completed checkpoint manifest has an unexpected resume offset for prefix: " + prefix);
        }
        if (have_segments && manifest.next_segment_id <= max_segment_id) {
            throw std::runtime_error("Checkpoint manifest next_segment_id does not advance past committed segments for prefix: " + prefix);
        }
    }

    void validate_segment_files(const string& genotype_path,
                                const string& predicted_path,
                                const string& variants_path,
                                const string& loo_path,
                                const vcf::PredictMissingCheckpointSegment& segment,
                                size_t sample_count) {
        if (!local_file_exists(genotype_path) || !local_file_exists(predicted_path) ||
            !local_file_exists(variants_path) || !local_file_exists(loo_path)) {
            throw std::runtime_error("referenced segment files are missing");
        }

        const uint64_t expected_genotype_bytes = static_cast<uint64_t>(segment.rows) * sample_count * sizeof(float);
        const uint64_t expected_predicted_bytes = static_cast<uint64_t>(segment.rows) * sample_count;
        if (local_file_size(genotype_path) != expected_genotype_bytes) {
            throw std::runtime_error("genotype checkpoint segment does not match manifest row count");
        }
        if (local_file_size(predicted_path) != expected_predicted_bytes) {
            throw std::runtime_error("predicted checkpoint segment does not match manifest row count");
        }
        if (count_lines(variants_path) != segment.rows) {
            throw std::runtime_error("variants checkpoint segment does not match manifest row count");
        }
        if (count_lines(loo_path) != segment.loo_rows) {
            throw std::runtime_error("loo checkpoint segment does not match manifest row count");
        }
    }

    string format_double(double value) {
        std::ostringstream out;
        out.precision(17);
        out << value;
        return out.str();
    }

    string serialize_loo_row(const vcf::ImputationLooRow& row) {
        std::ostringstream out;
        out << row.variant << '\t'
            << row.logical_index << '\t'
            << row.n_observed << '\t'
            << row.n_missing << '\t'
            << format_double(row.oob_mae) << '\t'
            << format_double(row.oob_rmse) << '\t'
            << format_double(row.rounded_acc);
        return out.str();
    }

    vcf::ImputationLooRow parse_loo_row(const string& line) {
        std::vector<string> tokens;
        std::istringstream input(line);
        string token;
        while (std::getline(input, token, '\t')) {
            tokens.push_back(token);
        }
        if (tokens.size() != 9) {
            throw std::runtime_error("Checkpoint loo segment row is malformed");
        }

        vcf::ImputationLooRow row;
        row.variant = tokens[0] + '\t' + tokens[1] + '\t' + tokens[2];
        row.logical_index = static_cast<size_t>(std::stoull(tokens[3]));
        row.n_observed = static_cast<size_t>(std::stoull(tokens[4]));
        row.n_missing = static_cast<size_t>(std::stoull(tokens[5]));
        row.oob_mae = std::stod(tokens[6]);
        row.oob_rmse = std::stod(tokens[7]);
        row.rounded_acc = std::stod(tokens[8]);
        return row;
    }
}

namespace vcf {
    PredictMissingCheckpoint::PredictMissingCheckpoint(std::string prefix)
        : prefix_value(std::move(prefix)),
          manifest_file(make_path(prefix_value, "_predict_state.tsv")),
          manifest_previous_file(make_path(prefix_value, "_predict_state.previous.tsv")),
          manifest_tmp_file(make_path(prefix_value, "_predict_state.tsv.tmp")) {}

    const string& PredictMissingCheckpoint::prefix() const {
        return prefix_value;
    }

    string PredictMissingCheckpoint::segment_genotype_file(size_t segment_id) const {
        return make_segment_path(prefix_value, segment_id, "_genotype.bin");
    }

    string PredictMissingCheckpoint::segment_predicted_file(size_t segment_id) const {
        return make_segment_path(prefix_value, segment_id, "_predicted.bin");
    }

    string PredictMissingCheckpoint::segment_variants_file(size_t segment_id) const {
        return make_segment_path(prefix_value, segment_id, "_variants.txt");
    }

    string PredictMissingCheckpoint::segment_loo_file(size_t segment_id) const {
        return make_segment_path(prefix_value, segment_id, "_loo.txt");
    }

    bool PredictMissingCheckpoint::has_manifest() const {
        return local_file_exists(manifest_file) || local_file_exists(manifest_previous_file);
    }

    void PredictMissingCheckpoint::reset() const {
        auto cleanup_manifest = [&](const string& path) {
            std::ifstream input(path);
            if (!input.good()) {
                return;
            }
            std::ostringstream text;
            text << input.rdbuf();
            try {
                auto manifest = parse_manifest_text(text.str());
                for (const auto& segment : manifest.segments) {
                    remove_if_exists(segment_genotype_file(segment.id));
                    remove_if_exists(segment_predicted_file(segment.id));
                    remove_if_exists(segment_variants_file(segment.id));
                    remove_if_exists(segment_loo_file(segment.id));
                }
            } catch (...) {
            }
        };

        cleanup_manifest(manifest_file);
        cleanup_manifest(manifest_previous_file);

        remove_if_exists(manifest_file);
        remove_if_exists(manifest_previous_file);
        remove_if_exists(manifest_tmp_file);
        sync_directory(manifest_file);
    }

    bool PredictMissingCheckpoint::load_manifest(PredictMissingCheckpointManifest& manifest,
                                                 std::size_t* visible_versions,
                                                 bool* recovered_from_previous,
                                                 std::string* diagnostics) const {
        const bool current_exists = local_file_exists(manifest_file);
        const bool previous_exists = local_file_exists(manifest_previous_file);
        if (visible_versions != nullptr) {
            *visible_versions = static_cast<size_t>(current_exists) + static_cast<size_t>(previous_exists);
        }
        if (recovered_from_previous != nullptr) {
            *recovered_from_previous = false;
        }
        if (diagnostics != nullptr) {
            diagnostics->clear();
        }

        auto load_candidate = [&](const string& path, const string& slot_name, PredictMissingCheckpointManifest& out) -> string {
            try {
                std::ifstream input(path);
                if (!input.good()) {
                    throw std::runtime_error("manifest file is not readable");
                }
                std::ostringstream buffer;
                buffer << input.rdbuf();
                out = parse_manifest_text(buffer.str());
                validate_manifest(out, prefix_value);
                for (const auto& segment : out.segments) {
                    validate_segment_files(segment_genotype_file(segment.id),
                                           segment_predicted_file(segment.id),
                                           segment_variants_file(segment.id),
                                           segment_loo_file(segment.id),
                                           segment,
                                           out.sample_count);
                }
                return string();
            } catch (const std::exception& e) {
                return slot_name + " checkpoint is unusable: " + e.what();
            }
        };

        string current_error;
        if (current_exists) {
            current_error = load_candidate(manifest_file, "current", manifest);
            if (current_error.empty()) {
                if (diagnostics != nullptr) {
                    *diagnostics = describe_manifest("current", manifest);
                }
                return true;
            }
        }

        PredictMissingCheckpointManifest previous_manifest;
        string previous_error;
        if (previous_exists) {
            previous_error = load_candidate(manifest_previous_file, "previous", previous_manifest);
            if (previous_error.empty()) {
                manifest = previous_manifest;
                if (recovered_from_previous != nullptr) {
                    *recovered_from_previous = true;
                }
                if (diagnostics != nullptr) {
                    *diagnostics = current_error.empty()
                        ? describe_manifest("previous", manifest)
                        : current_error + "; recovered using " + describe_manifest("previous", manifest);
                }
                return true;
            }
        }

        if (!current_exists && !previous_exists) {
            return false;
        }

        std::ostringstream error;
        error << "No usable checkpoint manifest found for prefix: " << prefix_value;
        if (!current_error.empty()) {
            error << "; " << current_error;
        }
        if (!previous_error.empty()) {
            error << "; " << previous_error;
        }
        throw std::runtime_error(error.str());
    }

    void PredictMissingCheckpoint::save_manifest(const PredictMissingCheckpointManifest& manifest) const {
        write_manifest_atomically(manifest_tmp_file, manifest_file, manifest_previous_file, manifest, true);
    }

    void PredictMissingCheckpoint::replace_current_manifest(const PredictMissingCheckpointManifest& manifest) const {
        write_manifest_atomically(manifest_tmp_file, manifest_file, manifest_previous_file, manifest, false);
    }

    PredictMissingCheckpointSegment PredictMissingCheckpoint::write_segment(
            std::size_t segment_id,
            const std::vector<Variant>& variants,
            const std::vector<std::vector<float>>& gmatrix,
            const std::vector<std::vector<bool>>& missing,
            const std::vector<ImputationLooRow>& loo_rows,
            std::size_t from_index,
            std::size_t to_index,
            std::size_t row_offset,
            std::size_t loo_from_index,
            std::size_t loo_to_index) const {
        if (to_index <= from_index && loo_to_index <= loo_from_index) {
            return {segment_id, 0, 0};
        }

        const string genotype_file = segment_genotype_file(segment_id);
        const string predicted_file = segment_predicted_file(segment_id);
        const string variants_file = segment_variants_file(segment_id);
        const string loo_file = segment_loo_file(segment_id);

        std::ofstream genotype_output(genotype_file, std::ios::binary | std::ios::trunc);
        std::ofstream predicted_output(predicted_file, std::ios::binary | std::ios::trunc);
        std::ofstream variants_output(variants_file, std::ios::trunc);
        std::ofstream loo_output(loo_file, std::ios::trunc);
        if (!genotype_output.good() || !predicted_output.good() || !variants_output.good() || !loo_output.good()) {
            throw std::runtime_error("Failed to create checkpoint segment for prefix: " + prefix_value);
        }

        std::vector<uint8_t> predicted_row;
        for (size_t absolute = from_index; absolute < to_index; ++absolute) {
            size_t local = absolute - row_offset;
            if (local >= gmatrix.size() || local >= missing.size() || local >= variants.size()) {
                throw std::runtime_error("Checkpoint write_segment received an out-of-range row index");
            }

            const auto& genotype_row = gmatrix[local];
            genotype_output.write(reinterpret_cast<const char*>(genotype_row.data()),
                                  static_cast<std::streamsize>(sizeof(float) * genotype_row.size()));

            predicted_row.assign(missing[local].begin(), missing[local].end());
            predicted_output.write(reinterpret_cast<const char*>(predicted_row.data()),
                                   static_cast<std::streamsize>(predicted_row.size()));

            variants_output << static_cast<std::string>(variants[local]) << "\n";
        }

        for (size_t index = loo_from_index; index < loo_to_index; ++index) {
            if (index >= loo_rows.size()) {
                throw std::runtime_error("Checkpoint write_segment received an out-of-range loo row index");
            }
            loo_output << serialize_loo_row(loo_rows[index]) << "\n";
        }

        genotype_output.flush();
        predicted_output.flush();
        variants_output.flush();
        loo_output.flush();
        if (!genotype_output.good() || !predicted_output.good() || !variants_output.good() || !loo_output.good()) {
            throw std::runtime_error("Failed to flush checkpoint segment for prefix: " + prefix_value);
        }
        genotype_output.close();
        predicted_output.close();
        variants_output.close();
        loo_output.close();
        if (!genotype_output.good() || !predicted_output.good() || !variants_output.good() || !loo_output.good()) {
            throw std::runtime_error("Failed to close checkpoint segment for prefix: " + prefix_value);
        }

        sync_file(genotype_file);
        sync_file(predicted_file);
        sync_file(variants_file);
        sync_file(loo_file);
        sync_directory(genotype_file);

        return {segment_id, to_index - from_index, loo_to_index - loo_from_index};
    }

    void PredictMissingCheckpoint::load_prefix(std::size_t n_samples,
                                               std::vector<std::string>& variant_names,
                                               std::vector<float>& genotype_values,
                                               std::vector<uint8_t>& predicted_values) const {
        PredictMissingCheckpointManifest manifest;
        if (!load_manifest(manifest) || manifest.flushed_rows == 0) {
            variant_names.clear();
            genotype_values.clear();
            predicted_values.clear();
            return;
        }
        if (manifest.sample_count != n_samples) {
            throw std::runtime_error("Checkpoint sample_count does not match the current request");
        }

        const size_t total_rows = manifest.flushed_rows;
        variant_names.clear();
        genotype_values.assign(total_rows * n_samples, 0.0f);
        predicted_values.assign(total_rows * n_samples, 0);

        size_t loaded_rows = 0;
        for (const auto& segment : manifest.segments) {
            std::ifstream variants_input(segment_variants_file(segment.id));
            std::ifstream genotype_input(segment_genotype_file(segment.id), std::ios::binary);
            std::ifstream predicted_input(segment_predicted_file(segment.id), std::ios::binary);
            if (!variants_input.good() || !genotype_input.good() || !predicted_input.good()) {
                throw std::runtime_error("Checkpoint segment files are incomplete for prefix: " + prefix_value);
            }

            const size_t offset = loaded_rows * n_samples;
            const size_t values = segment.rows * n_samples;
            genotype_input.read(reinterpret_cast<char*>(genotype_values.data() + offset),
                                static_cast<std::streamsize>(sizeof(float) * values));
            predicted_input.read(reinterpret_cast<char*>(predicted_values.data() + offset),
                                 static_cast<std::streamsize>(values));

            if (!genotype_input.good() || !predicted_input.good()) {
                throw std::runtime_error("Checkpoint binary segment does not match manifest row count");
            }

            string line;
            for (size_t row = 0; row < segment.rows; ++row) {
                if (!std::getline(variants_input, line)) {
                    throw std::runtime_error("Checkpoint variants segment does not match manifest row count");
                }
                variant_names.push_back(line);
            }
            loaded_rows += segment.rows;
        }

        if (loaded_rows != total_rows || variant_names.size() != total_rows) {
            throw std::runtime_error("Checkpoint prefix rows do not match the committed manifest");
        }
    }

    void PredictMissingCheckpoint::load_loo(std::vector<ImputationLooRow>& rows) const {
        PredictMissingCheckpointManifest manifest;
        if (!load_manifest(manifest) || manifest.flushed_loo_rows == 0) {
            rows.clear();
            return;
        }

        rows.clear();
        rows.reserve(manifest.flushed_loo_rows);
        for (const auto& segment : manifest.segments) {
            std::ifstream input(segment_loo_file(segment.id));
            if (!input.good()) {
                throw std::runtime_error("Checkpoint loo segment is missing for prefix: " + prefix_value);
            }
            string line;
            size_t segment_rows = 0;
            while (std::getline(input, line)) {
                rows.push_back(parse_loo_row(line));
                ++segment_rows;
            }
            if (segment_rows != segment.loo_rows) {
                throw std::runtime_error("Checkpoint loo segment does not match manifest row count");
            }
        }

        if (rows.size() != manifest.flushed_loo_rows) {
            throw std::runtime_error("Checkpoint loo rows do not match the committed manifest");
        }
    }
}