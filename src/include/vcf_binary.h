#ifndef SRC_VCF_BINARY_H
#define SRC_VCF_BINARY_H

#include <string>
#include <fstream>

#include "vcf_primitives.h"
#include "mio.hpp"

namespace vcf {
    class BinaryVCFScanner {
    public:
        virtual BinaryAllele scan(int pos) = 0;
    };

    class MemoryMappedScanner : public BinaryVCFScanner {
        mio::mmap_source mmap;
    public:
        explicit MemoryMappedScanner(const std::string& filename);
        BinaryAllele scan(int pos) override;
    };

    class IOScanner : public BinaryVCFScanner {
        std::ifstream in;
    public:
        explicit IOScanner(const std::string& filename);
        BinaryAllele scan(int pos) override;
    };
}

#endif //SRC_VCF_BINARY_H
