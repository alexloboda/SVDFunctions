#include "include/vcf_handlers.h"
#include <algorithm>
#include <iostream>

namespace {
    using std::vector;
}

namespace vcf {
    VariantsHandler::VariantsHandler(const std::vector<std::string>& samples) :samples(samples){}

    void VariantsHandler::processVariant(Variant variant, std::vector<Allele> alleles) {}

    bool VariantsHandler::isOfInterest(const Position& position) {
        return false;
    }

    CallRateHandler::CallRateHandler(const std::vector<std::string>& _samples, const std::vector<Range>& _ranges)
        :VariantsHandler(_samples), ranges(_ranges) {
        auto val = vector<int>();
        val.resize(samples.size());
        call_rate_matrix.resize(ranges.size(), val);
        n_variants.resize(ranges.size(), 0);
        std::sort(ranges.begin(), ranges.end());
    }

    void CallRateHandler::processVariant(Variant variant, std::vector<Allele> alleles) {
        for (int r = 0; r < ranges.size(); r++) {
            const Range& range = ranges[r];
            if (range.includes(variant.position())) {
                n_variants[r]++;
                for (int i = 0; i < alleles.size(); i++) {
                    if (alleles[i].alleleType() != MISSING) {
                        ++call_rate_matrix[r][i];
                    }
                }
            }
        }
    }

    bool CallRateHandler::isOfInterest(const Position& position) {
        auto it = std::lower_bound(ranges.begin(), ranges.end(), position);
        if (it == ranges.end()) {
            return false;
        }
        return it->includes(position);
    }

    void GenotypeMatrixHandler::processVariant(Variant variant, std::vector<Allele> alleles) {
        vector<AlleleType> row;
        for (const Allele& allele: alleles) {
            row.push_back(allele.alleleType());
        }
        gmatrix.push_back(row);
        variants.push_back(variant);
    }

    bool GenotypeMatrixHandler::isOfInterest(const Position& position) {
        return true;
    }

    BinaryFileHandler::BinaryFileHandler(const std::vector<std::string>& samples, std::string main_filename,
                                         std::string metadata_file) :VariantsHandler(samples),
                                         binary(main_filename, std::ios::binary), meta(metadata_file) {
        for (const std::string& sample: samples) {
            meta << sample << DELIM;
        }
        meta << "\n";
    }

    void BinaryFileHandler::processVariant(Variant variant, std::vector<Allele> alleles) {
        meta << (std::string)variant << "\n";
        for (const Allele& allele: alleles) {
            binary << BinaryAllele::fromAllele(allele);
        }
    }

    bool BinaryFileHandler::isOfInterest(const Position& position) {
        return true;
    }
}
