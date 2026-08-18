#ifndef SRC_CONTROLS_MATRIX_H
#define SRC_CONTROLS_MATRIX_H

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

// On-disk control matrices for the file backed control selection.
//
// Rows are variants and every row is stored contiguously, so a pass that only
// needs the case variants can seek past the rows it does not want, and a row
// can be decoded without materialising the rest of the matrix. That is the one
// property the whole file backed path depends on: both the projection and the
// per-cluster counting are single streaming passes over variant rows.
//
// Layout (little endian, as the package already assumes for the VCF binaries):
//
//   "SVDCTL"                        6 bytes
//   version                         uint32
//   dtype                           uint8
//   nrow (variants)                 int64
//   ncol (samples or clusters)      int64
//   data_offset                     int64
//   nrow row names, ncol col names  uint32 length + bytes each
//   ncol column weights             int32 each
//   padding                         up to an 8 byte boundary
//   nrow rows of ncol elements      at data_offset
//
// A column weight is the number of control samples the column stands for: one
// for a genotype matrix, the cluster size for collapsed cluster counts. The
// matching stage needs those sizes and a collapsed file has no samples left to
// count, so they have to travel with the data.
namespace ctlm {

enum class DType : std::uint8_t {
    F64 = 0,            // imputed genotypes, 8 bytes per entry
    F32 = 1,            // imputed genotypes, 4 bytes per entry
    U8 = 2,             // raw genotypes: 0, 1, 2, or 255 for missing
    CLUSTER_COUNTS = 3  // hom-ref/het/hom-alt byte triples, already collapsed
};

extern const char MAGIC[];
constexpr std::size_t MAGIC_SIZE = 6;
constexpr std::uint32_t VERSION = 1;
constexpr std::uint8_t U8_MISSING = 255;

std::size_t element_size(DType dtype);
// Whether rows of this type decode to genotype values (as opposed to counts).
bool is_numeric(DType dtype);
const char* dtype_name(DType dtype);

struct Header {
    DType dtype = DType::F64;
    std::int64_t nrow = 0;
    std::int64_t ncol = 0;
    std::int64_t data_offset = 0;
    std::vector<std::string> rownames;
    std::vector<std::string> colnames;
    std::vector<std::int32_t> colweights;

    std::size_t row_bytes() const;
};

// Fills in data_offset and writes everything up to (and including) the padding.
void write_header(std::ostream& out, Header& header);
Header read_header(const std::string& path);

// Opens a control matrix and validates that its payload is complete. Rows are
// read through RowReader, one of which is needed per thread.
class Reader {
    Header hdr;
    std::string file;
public:
    explicit Reader(const std::string& path);

    const Header& header() const { return hdr; }
    const std::string& path() const { return file; }
};

// Reads rows into a buffer of its own. Deliberately not a memory map: mapped
// pages of a large matrix accumulate in the resident set as they are touched,
// which would defeat the point of streaming. A buffer bounds what the process
// holds to a single row, and leaves the caching to the kernel.
class RowReader {
    const Header* hdr;
    std::ifstream in;
    std::vector<char> buffer;
    std::int64_t position = -1;

    void seek_and_read(std::int64_t offset, std::size_t bytes);
public:
    explicit RowReader(const Reader& reader);

    // Reads columns [from, to) of row `r` and returns a pointer to those bytes,
    // valid until the next call.
    const char* raw(std::int64_t r, std::int64_t from, std::int64_t to);
    const char* raw(std::int64_t r) { return raw(r, 0, hdr->ncol); }
    // Decodes columns [from, to) of row `r` into `out`, which must hold
    // `to - from` doubles. Missing entries decode to NaN; the return value is
    // the first column that was missing, or -1 when the range was complete.
    std::int64_t decode_row(std::int64_t r, std::int64_t from, std::int64_t to, double* out);
};

}

#endif //SRC_CONTROLS_MATRIX_H
