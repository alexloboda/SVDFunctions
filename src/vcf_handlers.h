#ifndef SRC_VCF_TABLE_H
#define SRC_VCF_TABLE_H

#include <vector>
#include <string>
#include "vcf_primitives.h"

namespace vcf {
    class VariantsHandler {
        const std::vector<std::string>& samples;

    public:
        VariantsHandler(const std::vector<std::string>& samples);
        virtual void processVariant(Variant variant, std::vector<AlleleType > alleles) = 0;
    };

    class CallRateHandler: public VariantsHandler {
        const std::vector<Range>& ranges;
        std::vector<std::vector<int>> call_rate_matrix;
    public:
        CallRateHandler(const std::vector<std::string>& samples, const std::vector<Range>& ranges);
        void processVariant(Variant variant, std::vector<AlleleType> alleles) override;
    };

    class GenotypeMatrixHandler: public VariantsHandler {
        std::vector<std::vector<int>> gmatrix;
        std::vector<Variant> variants;
    public:
        using VariantsHandler::VariantsHandler;
        void processVariant(Variant variant, std::vector<AlleleType> alleles);
    };
}

#endif //SRC_VCF_TABLE_H
