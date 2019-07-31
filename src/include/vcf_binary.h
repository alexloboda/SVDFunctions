#ifndef SRC_VCF_BINARY_H
#define SRC_VCF_BINARY_H

#include <string>
#include <fstream>

#include "vcf_primitives.h"
#include "third-party/mio.hpp"

namespace vcf {
    class BinaryVCFScanner {
    public:
        virtual BinaryAllele scan(size_t pos) const = 0;
    };

    class MemoryMappedScanner : public BinaryVCFScanner {
        std::string filename;
        mio::mmap_source mmap;
    public:
        explicit MemoryMappedScanner(const std::string& filename);
        MemoryMappedScanner(const MemoryMappedScanner& other);
        BinaryAllele scan(size_t pos) const override;
    };
}

#endif //SRC_VCF_BINARY_H
