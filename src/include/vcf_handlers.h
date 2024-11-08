#ifndef SRC_VCF_TABLE_H
#define SRC_VCF_TABLE_H

#include <vector>
#include <string>
#include <fstream>
#include <unordered_set>
#include "vcf_primitives.h"
#include "vcf_stats.h"

namespace vcf {
    class AlleleVector;

    class VariantsHandler {
    protected:
        const std::vector<std::string> samples;

    public:
        explicit VariantsHandler(const std::vector<std::string>& samples);
        virtual ~VariantsHandler() = default;
        virtual void processVariant(const Variant& variant, std::shared_ptr<AlleleVector>& alleles);
        virtual bool isOfInterest(const Variant& position);
    };

    class CallRateHandler: public VariantsHandler {
    protected:
        std::vector<Range> ranges;
        std::vector<int> n_variants;
        std::vector<std::vector<int>> call_rate_matrix;
    public:
        CallRateHandler(const std::vector<std::string>& samples, const std::vector<Range>& ranges);
        void processVariant(const Variant& variant, std::shared_ptr<AlleleVector>& alleles) override;
        bool isOfInterest(const Variant& position) override;
    };

    class GenotypeMatrixHandler;

    class GenotypeMatrixIterator {
        std::size_t pos;
        GenotypeMatrixHandler& gh;
    public:
        explicit GenotypeMatrixIterator(GenotypeMatrixHandler& gh);
        Variant operator*();
        void set(std::vector<float> genotypes);
        GenotypeMatrixIterator& operator++();
        bool dereferencable();
    };

    class GenotypeMatrixHandler: public VariantsHandler {
        const double EPS = 1e-8;
    protected:
        std::vector<std::vector<float>> gmatrix;
        std::vector<std::vector<bool>> missing;
        std::vector<Variant> variants;
        std::unordered_map<Variant, bool> available_variants;
        VCFFilterStats& stats;
        double missing_rate_threshold;
    public:
        GenotypeMatrixHandler(const std::vector<std::string>& samples, const std::vector<Variant>& variants,
                              VCFFilterStats& stats, double missing_rate_threshold);
        void processVariant(const Variant& variant, std::shared_ptr<AlleleVector>& alleles) override;
        bool isOfInterest(const Variant& position) override;
        std::vector<Variant> desired_variants();
        GenotypeMatrixIterator iterator();
        friend class GenotypeMatrixIterator;
    };

    class BinaryFileHandler: public VariantsHandler {
        const std::string DELIM = "\t";

        std::ofstream binary;
        std::ofstream meta;
    public:
        BinaryFileHandler(const std::vector<std::string>& samples, std::string main_filename,
                std::string metadata_file);
        void processVariant(const Variant& variant, std::shared_ptr<AlleleVector>& alleles) override;
        bool isOfInterest(const Variant& position) override;
    };
}

#endif //SRC_VCF_TABLE_H
