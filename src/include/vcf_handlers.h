#ifndef SRC_VCF_TABLE_H
#define SRC_VCF_TABLE_H

#include <vector>
#include <string>
#include <fstream>
#include <unordered_set>
#include "vcf_primitives.h"

namespace vcf {
    class VariantsHandler {
    protected:
        const std::vector<std::string> samples;

    public:
        explicit VariantsHandler(const std::vector<std::string>& samples);
        virtual void processVariant(Variant variant, std::vector<Allele> alleles);
        virtual bool isOfInterest(const Variant& position);
    };

    class CallRateHandler: public VariantsHandler {
    protected:
        std::vector<Range> ranges;
        std::vector<int> n_variants;
        std::vector<std::vector<int>> call_rate_matrix;
    public:
        CallRateHandler(const std::vector<std::string>& samples, const std::vector<Range>& ranges);
        void processVariant(Variant variant, std::vector<Allele> alleles) override;
        bool isOfInterest(const Variant& position) override;
    };

    class GenotypeMatrixHandler: public VariantsHandler {
        const double MISSING_RATE_THRESHOLD = 0.1;
    protected:
        std::vector<std::vector<AlleleType>> gmatrix;
        std::vector<Variant> variants;
        std::unordered_set<Variant> available_variants;
    public:
        GenotypeMatrixHandler(const std::vector<std::string>& samples, const std::vector<Variant>& variants);
        void processVariant(Variant variant, std::vector<Allele> alleles) override;
        bool isOfInterest(const Variant& position) override;
    };

    class BinaryFileHandler: public VariantsHandler {
        const std::string DELIM = "\t";

        std::ofstream binary;
        std::ofstream meta;
    public:
        BinaryFileHandler(const std::vector<std::string>& samples, std::string main_filename,
                std::string metadata_file);
        void processVariant(Variant variant, std::vector<Allele> alleles) override;
        bool isOfInterest(const Variant& position) override;
    };
}

#endif //SRC_VCF_TABLE_H
