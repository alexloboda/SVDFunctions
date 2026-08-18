#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <thread>
#include <vector>

#include "include/controls_matrix.h"
#include "include/third-party/cxxpool.h"

using namespace Rcpp;

namespace {

const double NA_DOUBLE = std::numeric_limits<double>::quiet_NaN();

ctlm::DType parse_dtype(const std::string& name) {
    if (name == "f64") return ctlm::DType::F64;
    if (name == "f32") return ctlm::DType::F32;
    if (name == "u8") return ctlm::DType::U8;
    if (name == "clusterCounts") return ctlm::DType::CLUSTER_COUNTS;
    Rcpp::stop("Unknown control matrix element type: " + name);
}

std::vector<std::string> to_strings(const CharacterVector& v) {
    std::vector<std::string> out;
    out.reserve(v.size());
    for (R_xlen_t i = 0; i < v.size(); i++) {
        out.push_back(Rcpp::as<std::string>(v[i]));
    }
    return out;
}

// R matrices are column major, the file is variant major, so writing means
// transposing. Doing that element by element would miss cache on every entry of
// a large matrix, so rows are staged in blocks and each block is transposed in
// small tiles that stay resident while they are filled.
constexpr std::int64_t TRANSPOSE_TILE = 64;
constexpr std::size_t STAGING_BUDGET_BYTES = 64u << 20;

double source_value(const double* src, std::size_t i) {
    return src[i];
}

double source_value(const int* src, std::size_t i) {
    return src[i] == NA_INTEGER ? NA_DOUBLE : static_cast<double>(src[i]);
}

void encode(char* dst, double value, ctlm::DType dtype, std::int64_t row, std::int64_t col) {
    switch (dtype) {
        case ctlm::DType::F64: {
            std::memcpy(dst, &value, sizeof(double));
            break;
        }
        case ctlm::DType::F32: {
            float narrowed = static_cast<float>(value);
            std::memcpy(dst, &narrowed, sizeof(float));
            break;
        }
        case ctlm::DType::U8: {
            std::uint8_t stored;
            if (std::isnan(value)) {
                stored = ctlm::U8_MISSING;
            } else if (value != std::floor(value) || value < 0.0 || value > 2.0) {
                throw std::runtime_error("dtype \"u8\" only stores genotypes 0, 1, 2 or NA, but entry [" +
                                         std::to_string(row + 1) + ", " + std::to_string(col + 1) +
                                         "] is " + std::to_string(value));
            } else {
                stored = static_cast<std::uint8_t>(value);
            }
            std::memcpy(dst, &stored, sizeof(std::uint8_t));
            break;
        }
        default:
            throw std::runtime_error("Genotype matrices cannot be written as cluster counts");
    }
}

template<typename T>
void write_matrix_rows(std::ostream& out, const T* src, std::int64_t nrow, std::int64_t ncol,
                       ctlm::DType dtype) {
    std::size_t elem = ctlm::element_size(dtype);
    std::size_t row_bytes = static_cast<std::size_t>(ncol) * elem;
    std::int64_t block_rows = row_bytes == 0 ? nrow
                                             : static_cast<std::int64_t>(STAGING_BUDGET_BYTES / row_bytes);
    block_rows = std::max<std::int64_t>(1, std::min(block_rows, nrow));

    std::vector<char> staging(static_cast<std::size_t>(block_rows) * row_bytes);

    for (std::int64_t row0 = 0; row0 < nrow; row0 += block_rows) {
        std::int64_t rows = std::min(block_rows, nrow - row0);
        Rcpp::checkUserInterrupt();

        for (std::int64_t jt = 0; jt < ncol; jt += TRANSPOSE_TILE) {
            std::int64_t j_end = std::min(ncol, jt + TRANSPOSE_TILE);
            for (std::int64_t it = 0; it < rows; it += TRANSPOSE_TILE) {
                std::int64_t i_end = std::min(rows, it + TRANSPOSE_TILE);
                for (std::int64_t j = jt; j < j_end; j++) {
                    std::size_t column_base = static_cast<std::size_t>(j) * nrow + row0;
                    for (std::int64_t i = it; i < i_end; i++) {
                        double value = source_value(src, column_base + i);
                        char* dst = staging.data() + static_cast<std::size_t>(i) * row_bytes +
                                    static_cast<std::size_t>(j) * elem;
                        encode(dst, value, dtype, row0 + i, j);
                    }
                }
            }
        }

        out.write(staging.data(), static_cast<std::streamsize>(static_cast<std::size_t>(rows) * row_bytes));
        if (!out) {
            throw std::runtime_error("Failed to write control matrix data");
        }
    }
}

}

// [[Rcpp::export]]
void write_controls_matrix_cpp(SEXP matrix, const std::string& path, const std::string& dtype,
                               const CharacterVector& rownames, const CharacterVector& colnames) {
    try {
        ctlm::Header header;
        header.dtype = parse_dtype(dtype);
        header.nrow = rownames.size();
        header.ncol = colnames.size();
        header.rownames = to_strings(rownames);
        header.colnames = to_strings(colnames);
        // Every column of a genotype matrix stands for exactly one sample.
        header.colweights.assign(static_cast<std::size_t>(header.ncol), 1);

        std::ofstream out(path, std::ios::binary | std::ios::trunc);
        if (!out) {
            Rcpp::stop("Cannot open control matrix file for writing: " + path);
        }
        ctlm::write_header(out, header);

        if (TYPEOF(matrix) == REALSXP) {
            write_matrix_rows(out, REAL(matrix), header.nrow, header.ncol, header.dtype);
        } else if (TYPEOF(matrix) == INTSXP) {
            write_matrix_rows(out, INTEGER(matrix), header.nrow, header.ncol, header.dtype);
        } else {
            Rcpp::stop("Control matrices must be numeric or integer");
        }
        out.close();
        if (!out) {
            Rcpp::stop("Failed to close control matrix file: " + path);
        }
    } catch (std::exception& e) {
        Rcpp::stop(e.what());
    }
}

// [[Rcpp::export]]
List controls_matrix_info_cpp(const std::string& path) {
    try {
        ctlm::Header header = ctlm::read_header(path);
        List ret;
        ret["dtype"] = std::string(ctlm::dtype_name(header.dtype));
        ret["nrow"] = static_cast<double>(header.nrow);
        ret["ncol"] = static_cast<double>(header.ncol);
        ret["rownames"] = CharacterVector(header.rownames.begin(), header.rownames.end());
        ret["colnames"] = CharacterVector(header.colnames.begin(), header.colnames.end());
        ret["colweights"] = IntegerVector(header.colweights.begin(), header.colweights.end());
        return ret;
    } catch (std::exception& e) {
        Rcpp::stop(e.what());
    }
}

// Streams the control matrix and projects it into the case PCA-like space.
//
// The projection is transition %*% (G - controlsMean), which decomposes into one
// contribution per variant, so the matrix never has to exist as a whole: each
// variant row is decoded, centered and folded into the k x n_samples result.
// Work is split over variants and every worker keeps its own accumulator, which
// is what keeps each worker's reads sequential -- the accumulators are cheap
// precisely because the reduced space is small.
// [[Rcpp::export]]
NumericMatrix project_controls_file_cpp(const std::string& path, const IntegerVector& variant_rows,
                                        const NumericMatrix& transition,
                                        const NumericVector& controls_mean, int threads) {
    using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    constexpr std::size_t TILE_BUDGET_BYTES = 8u << 20;
    constexpr std::size_t ACCUMULATOR_BUDGET_BYTES = 512u << 20;
    constexpr std::int64_t MAX_SLOTS = 64;

    try {
        ctlm::Reader reader(path);
        const ctlm::Header& header = reader.header();
        if (!ctlm::is_numeric(header.dtype)) {
            Rcpp::stop("Projection needs a genotype matrix, but " + path + " holds cluster counts");
        }

        std::int64_t n_samples = header.ncol;
        std::int64_t n_variants = variant_rows.size();
        std::int64_t k = transition.nrow();
        if (transition.ncol() != n_variants) {
            Rcpp::stop("transition must have one column per case variant");
        }
        if (controls_mean.size() != n_variants) {
            Rcpp::stop("controlsMean must have one value per case variant");
        }
        if (n_samples == 0 || n_variants == 0 || k == 0) {
            return NumericMatrix(k, n_samples);
        }

        Eigen::Map<const Eigen::MatrixXd> transition_map(transition.begin(), k, n_variants);

        // The variants are split into reduction slots whose number depends only on
        // the problem size, never on how many threads are available. Floating point
        // addition is not associative, so a thread dependent grouping would make the
        // reduced coordinates depend on the thread count, and through them the
        // discrete search downstream -- exactly the guarantee selectControls makes
        // about its other threading knobs.
        std::size_t accumulator_bytes = static_cast<std::size_t>(k) * n_samples * sizeof(double);
        std::int64_t slots = accumulator_bytes == 0
                ? 1 : static_cast<std::int64_t>(ACCUMULATOR_BUDGET_BYTES / accumulator_bytes);
        slots = std::max<std::int64_t>(1, std::min<std::int64_t>({slots, MAX_SLOTS, n_variants}));
        std::int64_t per_slot = (n_variants + slots - 1) / slots;

        int pool_size = threads > 0 ? threads : static_cast<int>(std::thread::hardware_concurrency());
        pool_size = std::max(1, std::min<int>(pool_size, static_cast<int>(slots)));

        std::size_t row_bytes = static_cast<std::size_t>(n_samples) * sizeof(double);
        std::int64_t tile = row_bytes == 0 ? 1 : static_cast<std::int64_t>(TILE_BUDGET_BYTES / row_bytes);
        tile = std::max<std::int64_t>(1, std::min<std::int64_t>(tile, 64));

        auto run = [&](std::int64_t begin, std::int64_t end) -> Eigen::MatrixXd {
            ctlm::RowReader rows(reader);
            Eigen::MatrixXd accumulator = Eigen::MatrixXd::Zero(k, n_samples);
            RowMajorMatrix centered(tile, n_samples);
            for (std::int64_t t0 = begin; t0 < end; t0 += tile) {
                std::int64_t len = std::min(tile, end - t0);
                for (std::int64_t c = 0; c < len; c++) {
                    std::int64_t variant = t0 + c;
                    std::int64_t missing = rows.decode_row(variant_rows[variant], 0, n_samples,
                                                           centered.row(c).data());
                    if (missing >= 0) {
                        throw std::runtime_error("Control genotypes must not be missing, but case variant " +
                                                 std::to_string(variant + 1) + " is missing for sample " +
                                                 std::to_string(missing + 1));
                    }
                    centered.row(c).array() -= controls_mean[variant];
                }
                accumulator.noalias() += transition_map.middleCols(t0, len) * centered.topRows(len);
            }
            return accumulator;
        };

        Eigen::MatrixXd points = Eigen::MatrixXd::Zero(k, n_samples);
        if (pool_size == 1) {
            for (std::int64_t begin = 0; begin < n_variants; begin += per_slot) {
                points += run(begin, std::min(n_variants, begin + per_slot));
            }
        } else {
            cxxpool::thread_pool pool(pool_size);
            std::vector<std::future<Eigen::MatrixXd>> futures;
            for (std::int64_t begin = 0; begin < n_variants; begin += per_slot) {
                std::int64_t end = std::min(n_variants, begin + per_slot);
                futures.push_back(pool.push([&run, begin, end]() { return run(begin, end); }));
            }
            // Summed in slot order, which is what makes the total independent of
            // the order the workers happened to finish in.
            for (auto& future: futures) {
                points += future.get();
            }
        }

        NumericMatrix ret(k, n_samples);
        Eigen::Map<Eigen::MatrixXd>(ret.begin(), k, n_samples) = points;
        return ret;
    } catch (std::exception& e) {
        Rcpp::stop(e.what());
    }
}

// Collapses a raw genotype matrix into the per-cluster allele counts the matching
// stage consumes, one variant row at a time, so neither the source matrix nor the
// result has to be held in R.
// [[Rcpp::export]]
List collapse_controls_matrix_cpp(const std::string& path, const IntegerVector& cluster_ids,
                                  const CharacterVector& cluster_labels, const std::string& out_path) {
    try {
        ctlm::Reader reader(path);
        const ctlm::Header& header = reader.header();
        if (header.dtype != ctlm::DType::U8) {
            Rcpp::stop("Collapsing needs a raw genotype matrix written with dtype \"u8\"");
        }
        std::int64_t n_samples = header.ncol;
        if (cluster_ids.size() != n_samples) {
            Rcpp::stop("clusterIds must have one entry per column of the control matrix");
        }
        std::size_t n_clusters = cluster_labels.size();
        if (n_clusters == 0) {
            Rcpp::stop("At least one cluster is required");
        }

        std::vector<int> ids(cluster_ids.begin(), cluster_ids.end());
        std::vector<int> sizes(n_clusters, 0);
        for (int id: ids) {
            if (id < 0 || static_cast<std::size_t>(id) >= n_clusters) {
                Rcpp::stop("Cluster index out of range");
            }
            ++sizes[id];
        }
        for (std::size_t c = 0; c < n_clusters; c++) {
            if (sizes[c] > 255) {
                Rcpp::stop("Each control cluster must contain at most 255 samples");
            }
        }

        ctlm::Header out_header;
        out_header.dtype = ctlm::DType::CLUSTER_COUNTS;
        out_header.nrow = header.nrow;
        out_header.ncol = static_cast<std::int64_t>(n_clusters);
        out_header.rownames = header.rownames;
        out_header.colnames = to_strings(cluster_labels);
        out_header.colweights.assign(sizes.begin(), sizes.end());

        std::ofstream out(out_path, std::ios::binary | std::ios::trunc);
        if (!out) {
            Rcpp::stop("Cannot open cluster counts file for writing: " + out_path);
        }
        ctlm::write_header(out, out_header);

        ctlm::RowReader rows(reader);
        std::vector<std::uint8_t> block(n_clusters * 3);
        for (std::int64_t r = 0; r < header.nrow; r++) {
            if (r % 1000 == 0) {
                Rcpp::checkUserInterrupt();
            }
            std::fill(block.begin(), block.end(), 0);
            const auto* row = reinterpret_cast<const std::uint8_t*>(rows.raw(r));
            for (std::int64_t j = 0; j < n_samples; j++) {
                std::uint8_t value = row[j];
                if (value <= 2) {
                    block[static_cast<std::size_t>(ids[j]) * 3 + value] += 1;
                }
            }
            out.write(reinterpret_cast<const char*>(block.data()),
                      static_cast<std::streamsize>(block.size()));
            if (!out) {
                Rcpp::stop("Failed to write cluster counts data");
            }
        }
        out.close();

        List ret;
        ret["path"] = out_path;
        ret["variants"] = static_cast<double>(header.nrow);
        ret["clusters"] = static_cast<int>(n_clusters);
        ret["clusterSizes"] = IntegerVector(sizes.begin(), sizes.end());
        return ret;
    } catch (std::exception& e) {
        Rcpp::stop(e.what());
    }
}
