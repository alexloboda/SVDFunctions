#include "include/vcf_binary.h"

#include <fstream>

vcf::MemoryMappedScanner::MemoryMappedScanner(const std::string& filename) :filename(filename) {
    std::error_code error;
    mmap.map(filename, error);
    if (error) {
        throw ParserException("Cannot map given file");
    }
}

vcf::BinaryAllele vcf::MemoryMappedScanner::scan(size_t pos) const {
    std::unique_ptr<char[]> bytes{new char[sizeof(BinaryAllele)]};
    size_t offset = pos * sizeof(BinaryAllele);
    for (size_t i = 0; i < sizeof(BinaryAllele); i++) {
        bytes[i] = mmap[offset + i];
    }
    return *(BinaryAllele*)bytes.get();
}

vcf::MemoryMappedScanner::MemoryMappedScanner(const vcf::MemoryMappedScanner &other)
                                :MemoryMappedScanner(other.filename) {}
