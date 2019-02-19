#include "vcf_handlers.h"

namespace {
    using std::vector;
}

namespace vcf {
    VariantsHandler::VariantsHandler(const std::vector<std::string>& samples) :samples(samples){}

    CallRateHandler::CallRateHandler(const std::vector<std::string>& samples, const std::vector<Range>& ranges)
        :VariantsHandler(samples), ranges(ranges) {
        auto val = vector<int>();
        val.resize(samples.size());
        call_rate_matrix.resize(ranges.size(), val);
    }

    void CallRateHandler::processVariant(Variant variant, std::vector<AlleleType> alleles) {
        for (int r = 0; r < ranges.size(); r++) {
            const Range& range = ranges[r];
            if (range.includes(variant.position())) {
                for (int i = 0; i < alleles.size(); i++) {
                    if (alleles[i].alleleType() != MISSING) {
                        ++call_rate_matrix[r][i];
                    }
                }
            }
        }
    }

    void GenotypeMatrixHandler::processVariant(Variant variant, std::vector<AlleleType> alleles) {
        vector<int> row;
        for (const AlleleType& allele: alleles) {
            row.push_back((int)allele);
        }
        gmatrix.push_back(row);
        variants.push_back(variant);
    }
}
