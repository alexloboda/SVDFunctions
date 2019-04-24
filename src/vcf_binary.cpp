#include "include/vcf_binary.h"

#include <fstream>

vcf::MemoryMappedScanner::MemoryMappedScanner(const std::string& filename) {
    std::error_code error;
    mmap.map(filename, error);
    if (error) {
        throw ParserException("Cannot map given file");
    }
}

vcf::BinaryAllele vcf::MemoryMappedScanner::scan(int pos) {
    char* bytes = new char[sizeof(BinaryAllele)];
    int offset = pos * sizeof(BinaryAllele);
    for (size_t i = 0; i < sizeof(BinaryAllele); i++) {
        bytes[i] = mmap[offset + i];
    }
    auto* obj = (BinaryAllele*)bytes;
    return *obj;
}

vcf::IOScanner::IOScanner(const std::string& filename) :in(filename, std::ios_base::binary) {}

vcf::BinaryAllele vcf::IOScanner::scan(int pos) {
    in.seekg(pos * sizeof(BinaryAllele));
    BinaryAllele allele{};
    in >> allele;
    return allele;
}
