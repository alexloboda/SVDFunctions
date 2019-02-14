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
        virtual void processVariant(Variant variant, std::vector<Allele> alleles) = 0;
    };
}

#endif //SRC_VCF_TABLE_H
