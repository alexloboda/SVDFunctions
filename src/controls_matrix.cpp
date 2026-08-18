#include "include/controls_matrix.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <ostream>
#include <stdexcept>

namespace ctlm {

const char MAGIC[] = "SVDCTL";

namespace {

constexpr std::size_t FIXED_HEADER_BYTES =
        MAGIC_SIZE + sizeof(std::uint32_t) + sizeof(std::uint8_t) + 3 * sizeof(std::int64_t);

std::int64_t align8(std::int64_t offset) {
    return (offset + 7) & ~static_cast<std::int64_t>(7);
}

template<typename T>
void put(std::ostream& out, T value) {
    out.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

template<typename T>
T get(std::istream& in) {
    T value;
    in.read(reinterpret_cast<char*>(&value), sizeof(T));
    if (in.gcount() != static_cast<std::streamsize>(sizeof(T))) {
        throw std::runtime_error("Control matrix file is truncated");
    }
    return value;
}

std::int64_t names_bytes(const std::vector<std::string>& names) {
    std::int64_t total = 0;
    for (const std::string& name: names) {
        total += static_cast<std::int64_t>(sizeof(std::uint32_t) + name.size());
    }
    return total;
}

void write_names(std::ostream& out, const std::vector<std::string>& names) {
    for (const std::string& name: names) {
        put<std::uint32_t>(out, static_cast<std::uint32_t>(name.size()));
        out.write(name.data(), static_cast<std::streamsize>(name.size()));
    }
}

void write_weights(std::ostream& out, const std::vector<std::int32_t>& weights) {
    for (std::int32_t weight: weights) {
        put<std::int32_t>(out, weight);
    }
}

std::vector<std::int32_t> read_weights(std::istream& in, std::int64_t n) {
    std::vector<std::int32_t> weights;
    weights.reserve(static_cast<std::size_t>(n));
    for (std::int64_t i = 0; i < n; i++) {
        weights.push_back(get<std::int32_t>(in));
    }
    return weights;
}

std::vector<std::string> read_names(std::istream& in, std::int64_t n) {
    std::vector<std::string> names;
    names.reserve(static_cast<std::size_t>(n));
    for (std::int64_t i = 0; i < n; i++) {
        auto length = get<std::uint32_t>(in);
        std::string name(length, '\0');
        if (length > 0) {
            in.read(&name[0], length);
            if (in.gcount() != static_cast<std::streamsize>(length)) {
                throw std::runtime_error("Control matrix file is truncated");
            }
        }
        names.push_back(std::move(name));
    }
    return names;
}

}

std::size_t element_size(DType dtype) {
    switch (dtype) {
        case DType::F64: return 8;
        case DType::F32: return 4;
        case DType::U8: return 1;
        case DType::CLUSTER_COUNTS: return 3;
    }
    throw std::invalid_argument("Unknown control matrix element type");
}

bool is_numeric(DType dtype) {
    return dtype == DType::F64 || dtype == DType::F32 || dtype == DType::U8;
}

const char* dtype_name(DType dtype) {
    switch (dtype) {
        case DType::F64: return "f64";
        case DType::F32: return "f32";
        case DType::U8: return "u8";
        case DType::CLUSTER_COUNTS: return "clusterCounts";
    }
    throw std::invalid_argument("Unknown control matrix element type");
}

std::size_t Header::row_bytes() const {
    return static_cast<std::size_t>(ncol) * element_size(dtype);
}

void write_header(std::ostream& out, Header& header) {
    if (static_cast<std::int64_t>(header.colweights.size()) != header.ncol) {
        throw std::runtime_error("Control matrix needs one column weight per column");
    }
    std::int64_t offset = static_cast<std::int64_t>(FIXED_HEADER_BYTES) +
                          names_bytes(header.rownames) + names_bytes(header.colnames) +
                          header.ncol * static_cast<std::int64_t>(sizeof(std::int32_t));
    header.data_offset = align8(offset);

    out.write(MAGIC, MAGIC_SIZE);
    put<std::uint32_t>(out, VERSION);
    put<std::uint8_t>(out, static_cast<std::uint8_t>(header.dtype));
    put<std::int64_t>(out, header.nrow);
    put<std::int64_t>(out, header.ncol);
    put<std::int64_t>(out, header.data_offset);
    write_names(out, header.rownames);
    write_names(out, header.colnames);
    write_weights(out, header.colweights);

    static const char padding[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    out.write(padding, header.data_offset - offset);
}

Header read_header(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Cannot open control matrix file: " + path);
    }

    char magic[MAGIC_SIZE];
    in.read(magic, MAGIC_SIZE);
    if (in.gcount() != static_cast<std::streamsize>(MAGIC_SIZE) ||
            std::memcmp(magic, MAGIC, MAGIC_SIZE) != 0) {
        throw std::runtime_error(path + " is not a control matrix file");
    }

    auto version = get<std::uint32_t>(in);
    if (version != VERSION) {
        throw std::runtime_error("Unsupported control matrix version in " + path);
    }

    Header header;
    auto dtype = get<std::uint8_t>(in);
    if (dtype > static_cast<std::uint8_t>(DType::CLUSTER_COUNTS)) {
        throw std::runtime_error("Unknown control matrix element type in " + path);
    }
    header.dtype = static_cast<DType>(dtype);
    header.nrow = get<std::int64_t>(in);
    header.ncol = get<std::int64_t>(in);
    header.data_offset = get<std::int64_t>(in);
    if (header.nrow < 0 || header.ncol < 0) {
        throw std::runtime_error("Negative dimensions in " + path);
    }
    header.rownames = read_names(in, header.nrow);
    header.colnames = read_names(in, header.ncol);
    header.colweights = read_weights(in, header.ncol);
    return header;
}

Reader::Reader(const std::string& path) : hdr(read_header(path)), file(path) {
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    if (!in) {
        throw std::runtime_error("Cannot open control matrix file: " + path);
    }
    auto size = static_cast<std::uint64_t>(in.tellg());
    std::uint64_t expected = static_cast<std::uint64_t>(hdr.data_offset) +
                             static_cast<std::uint64_t>(hdr.nrow) * hdr.row_bytes();
    if (size < expected) {
        throw std::runtime_error("Control matrix file is truncated: " + path);
    }
}

RowReader::RowReader(const Reader& reader)
        : hdr(&reader.header()), in(reader.path(), std::ios::binary),
          buffer(hdr->row_bytes()) {
    if (!in) {
        throw std::runtime_error("Cannot open control matrix file: " + reader.path());
    }
}

void RowReader::seek_and_read(std::int64_t offset, std::size_t bytes) {
    // Rows are usually consumed in order, so skip the seek when the handle is
    // already where it needs to be.
    if (offset != position) {
        in.seekg(offset, std::ios::beg);
        if (!in) {
            throw std::runtime_error("Cannot seek in control matrix file");
        }
    }
    in.read(buffer.data(), static_cast<std::streamsize>(bytes));
    if (in.gcount() != static_cast<std::streamsize>(bytes)) {
        throw std::runtime_error("Control matrix file is truncated");
    }
    position = offset + static_cast<std::int64_t>(bytes);
}

const char* RowReader::raw(std::int64_t r, std::int64_t from, std::int64_t to) {
    if (r < 0 || r >= hdr->nrow) {
        throw std::out_of_range("Control matrix row out of range");
    }
    if (from < 0 || to > hdr->ncol || from > to) {
        throw std::out_of_range("Control matrix column range out of bounds");
    }
    std::size_t elem = element_size(hdr->dtype);
    std::int64_t offset = hdr->data_offset +
                          static_cast<std::int64_t>(r) * static_cast<std::int64_t>(hdr->row_bytes()) +
                          from * static_cast<std::int64_t>(elem);
    seek_and_read(offset, static_cast<std::size_t>(to - from) * elem);
    return buffer.data();
}

std::int64_t RowReader::decode_row(std::int64_t r, std::int64_t from, std::int64_t to, double* out) {
    if (!is_numeric(hdr->dtype)) {
        throw std::invalid_argument("Control matrix does not hold genotypes");
    }

    const char* base = raw(r, from, to);
    std::int64_t missing = -1;
    const double NaN = std::numeric_limits<double>::quiet_NaN();

    switch (hdr->dtype) {
        case DType::F64: {
            const char* p = base;
            for (std::int64_t j = from; j < to; j++, p += sizeof(double)) {
                double value;
                std::memcpy(&value, p, sizeof(double));
                out[j - from] = value;
                if (missing < 0 && std::isnan(value)) {
                    missing = j;
                }
            }
            break;
        }
        case DType::F32: {
            const char* p = base;
            for (std::int64_t j = from; j < to; j++, p += sizeof(float)) {
                float value;
                std::memcpy(&value, p, sizeof(float));
                out[j - from] = value;
                if (missing < 0 && std::isnan(value)) {
                    missing = j;
                }
            }
            break;
        }
        default: {
            const auto* p = reinterpret_cast<const std::uint8_t*>(base);
            for (std::int64_t j = from; j < to; j++, p++) {
                if (*p == U8_MISSING) {
                    out[j - from] = NaN;
                    if (missing < 0) {
                        missing = j;
                    }
                } else {
                    out[j - from] = *p;
                }
            }
            break;
        }
    }

    return missing;
}

}
