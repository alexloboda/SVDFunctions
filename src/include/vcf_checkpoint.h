#ifndef SRC_VCF_CHECKPOINT_H
#define SRC_VCF_CHECKPOINT_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "vcf_primitives.h"

namespace vcf {
    struct ImputationLooRow;

    struct PredictMissingCheckpointSegment {
        std::size_t id = 0;
        std::size_t rows = 0;
        std::size_t loo_rows = 0;
    };

    struct PredictMissingCheckpointManifest {
        std::size_t generation = 0;
        uint64_t signature = 0;
        std::size_t sample_count = 0;
        std::size_t flushed_rows = 0;
        std::size_t flushed_loo_rows = 0;
        int64_t resume_offset = -1;
        std::size_t resume_row_index = 0;
        bool complete = false;
        std::size_t next_segment_id = 0;
        std::vector<PredictMissingCheckpointSegment> segments;
    };

    class PredictMissingCheckpoint {
        std::string prefix_value;
        std::string manifest_file;
        std::string manifest_previous_file;
        std::string manifest_tmp_file;

        std::string segment_genotype_file(std::size_t segment_id) const;
        std::string segment_predicted_file(std::size_t segment_id) const;
        std::string segment_variants_file(std::size_t segment_id) const;
        std::string segment_loo_file(std::size_t segment_id) const;
    public:
        explicit PredictMissingCheckpoint(std::string prefix);

        const std::string& prefix() const;
        bool has_manifest() const;
        void reset() const;
        bool load_manifest(PredictMissingCheckpointManifest& manifest,
                           std::size_t* visible_versions = nullptr,
                           bool* recovered_from_previous = nullptr,
                           std::string* diagnostics = nullptr) const;
        void save_manifest(const PredictMissingCheckpointManifest& manifest) const;
        void replace_current_manifest(const PredictMissingCheckpointManifest& manifest) const;
        PredictMissingCheckpointSegment write_segment(std::size_t segment_id,
                                                     const std::vector<Variant>& variants,
                                                     const std::vector<std::vector<float>>& gmatrix,
                                                     const std::vector<std::vector<bool>>& missing,
                                                     const std::vector<ImputationLooRow>& loo_rows,
                                                     std::size_t from_index,
                                                     std::size_t to_index,
                                                     std::size_t row_offset,
                                                     std::size_t loo_from_index,
                                                     std::size_t loo_to_index) const;
        void load_prefix(std::size_t n_samples,
                         std::vector<std::string>& variant_names,
                         std::vector<float>& genotype_values,
                         std::vector<uint8_t>& predicted_values) const;
        void load_loo(std::vector<ImputationLooRow>& rows) const;
    };
}

#endif //SRC_VCF_CHECKPOINT_H