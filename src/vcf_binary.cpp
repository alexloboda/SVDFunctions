#include "include/vcf_binary.h"

#include <fstream>

vcf::MemoryMappedScanner::MemoryMappedScanner(const std::string& filename) {
    std::error_code error;
    mmap.map(filename, error);
    if (error) {
        throw ParserException("Cannot map given file");
    }
}

vcf::BinaryAllele vcf::MemoryMappedScanner::scan(size_t pos) const {
    char* bytes = new char[sizeof(BinaryAllele)];
    size_t offset = pos * sizeof(BinaryAllele);
    for (size_t i = 0; i < sizeof(BinaryAllele); i++) {
        bytes[i] = mmap[offset + i];
    }
    auto* obj = (BinaryAllele*)bytes;
    return *obj;
}
