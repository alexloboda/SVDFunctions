#include "include/vcf_handlers.h"
#include <algorithm>
#include <iostream>

namespace {
    using std::vector;
    using std::any_of;

    bool MAC_filter(const std::vector<vcf::Allele>& alleles) {
        auto ref_f = [](const vcf::Allele& a) {return a.alleleType() == vcf::HOMREF || a.alleleType() == vcf::HET;};
        auto alt_f = [](const vcf::Allele& a) {return a.alleleType() == vcf::HOM || a.alleleType() == vcf::HET;};
        return any_of(alleles.begin(), alleles.end(), ref_f) && any_of(alleles.begin(), alleles.end(), alt_f);
    }
}

namespace vcf {
    VariantsHandler::VariantsHandler(const std::vector<std::string>& samples) :samples(samples){}

    void VariantsHandler::processVariant(const Variant& variant, const std::vector<Allele>& alleles) {}

    bool VariantsHandler::isOfInterest(const Variant& variant) {
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

    void CallRateHandler::processVariant(const Variant& variant, const std::vector<Allele>& alleles) {
        auto it = std::lower_bound(ranges.begin(), ranges.end(), variant.position());
        if (it == ranges.end() || !it->includes(variant.position())) {
            return;
        }
        auto r = std::distance(ranges.begin(), it);

        n_variants[r]++;
        for (int i = 0; i < alleles.size(); i++) {
            if (alleles[i].alleleType() != MISSING) {
                ++call_rate_matrix[r][i];
            }
        }
    }

    bool CallRateHandler::isOfInterest(const Variant& variant) {
        const Position& position = variant.position();
        auto it = std::lower_bound(ranges.begin(), ranges.end(), position);
        if (it == ranges.end()) {
            return false;
        }
        return it->includes(position);
    }

    void GenotypeMatrixHandler::processVariant(const Variant& variant, const std::vector<Allele>& alleles) {
        if (!isOfInterest(variant)) {
            return;
        }
        double missing_rate = std::count_if(alleles.begin(), alleles.end(), [](const Allele& x){
            return x.alleleType() == MISSING;
        }) / (double)alleles.size();
        if (missing_rate > MISSING_RATE_THRESHOLD) {
            stats.add(Stat::MISSING_RATE, 1);
            return;
        }
        if (!MAC_filter(alleles)) {
            return;
        }
        vector<AlleleType> row;
        for (const Allele& allele: alleles) {
            row.push_back(allele.alleleType());
        }
        gmatrix.push_back(row);
        variants.push_back(variant);
    }

    bool GenotypeMatrixHandler::isOfInterest(const Variant& variant) {
        return available_variants.empty() || available_variants.find(variant) != available_variants.end();
    }

    GenotypeMatrixHandler::GenotypeMatrixHandler(const std::vector<std::string>& ss, const std::vector<Variant>& vs,
                                                 VCFFilterStats& stats)
        : VariantsHandler(ss), stats(stats) {
        available_variants.insert(vs.begin(), vs.end());
    }

    BinaryFileHandler::BinaryFileHandler(const std::vector<std::string>& samples, std::string main_filename,
                                         std::string metadata_file) :VariantsHandler(samples),
                                         binary(main_filename, std::ios::binary), meta(metadata_file) {
        for (const std::string& sample: samples) {
            meta << sample << DELIM;
        }
        meta << "\n";
    }

    void BinaryFileHandler::processVariant(const Variant& variant, const std::vector<Allele>& alleles) {
        if (!MAC_filter(alleles)) {
            return;
        }
        meta << (std::string)variant << "\n";
        for (const Allele& allele: alleles) {
            binary << BinaryAllele::fromAllele(allele);
        }
    }

    bool BinaryFileHandler::isOfInterest(const Variant& variant) {
        return true;
    }
}
