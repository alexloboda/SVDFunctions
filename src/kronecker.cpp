#include "include/kronecker.h"

#include <fstream>
#include "include/third-party/cxxpool.h"
#include "include/third-party/irlba/irlba.hpp"

#include "include/third-party/zstr/zstr.hpp"

namespace {
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
}

namespace {
    
struct nkp_result {
    Matrix B;
    Matrix C;
};

nkp_result nkp(const Matrix& A, int n1, int m1, int n2, int m2) {
    int n = A.rows();
    int m = A.cols();
    
    assert(n1 == m1);
    assert(n2 == m2);

    Matrix A_mod = 0.5 * (A + A.transpose());

    Matrix R(n1 * m1, n2 * m2);
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < m1; ++j) {
            for (int k = 0; k < n2; ++k) {
                for (int l = 0; l < m2; ++l) {
                    R(i * m1 + j, k * m2 + l) = A_mod(k * n1 + i, l * m1 + j);
                }
            }
        }
    }

    // SVD

    irlba::Options opt;
    opt.extra_work = 20;
    opt.max_iterations = 50;
    auto svd = irlba::compute(R, 1, opt);;

    Matrix B = svd.U;
    double S = svd.D(0, 0);
    Matrix C = svd.V;

    double SqrtS = std::sqrt(S);

    B = B * SqrtS;
    C = C * SqrtS;

    B.resize(n1, m1);
    C.resize(n2, m2);

    B = 0.5 * (B + B.transpose());
    C = 0.5 * (C + C.transpose());

    if ((B.diagonal().array() < 0).all() && (C.diagonal().array() < 0).all()) {
        B = -B;
        C = -C;
    }

    return {B, C};
}

std::vector<Matrix> preprocess(const Matrix& lhs, const Matrix& rhs, int n) {
    Matrix ds(lhs.rows() * rhs.rows(), lhs.cols());
    for (int i = 0; i < lhs.rows(); i++) {
        for (int j = 0; j < rhs.rows(); j++) {
            ds.row(i * rhs.rows() + j) = lhs.row(i) - rhs.row(j);
        }
    }

    std::vector<Matrix> ret;

    Matrix prod = ds;
    for (int i = 0; i < n - 1; i++) {
        prod = Eigen::kroneckerProduct(prod, ds).eval();
    }


}
}

namespace {
using Eigen::placeholders::all;

template<class Derived>
void write_binary(std::ostream &os, const Eigen::PlainObjectBase<Derived> &matrix)
{
    typedef typename Derived::Index Index;
    typedef typename Derived::Scalar Scalar;

    assert(matrix.rows() == matrix.cols() && "Matrix must be square for symmetry.");

    Index size = matrix.rows();
    os.write((char*) (&size), sizeof(Index));

    // Write only the upper triangular part
    for (Index i = 0; i < size; ++i)
        for (Index j = i; j < size; ++j)
            os.write((char*) (&matrix(i, j)), sizeof(Scalar));
}

template<class Derived>
void read_binary(std::istream &is, Eigen::PlainObjectBase<Derived> &matrix)
{
    typedef typename Derived::Index Index;
    typedef typename Derived::Scalar Scalar;

    Index size = 0;
    is.read((char*) (&size), sizeof(Index));

    matrix.resize(size, size);

    // Read and reconstruct the matrix from the upper triangular part
    for (Index i = 0; i < size; ++i)
        for (Index j = i; j < size; ++j)
        {
            Scalar value;
            is.read((char*) (&value), sizeof(Scalar));
            matrix(i, j) = value;
            if (i != j)
                matrix(j, i) = value; // Fill the symmetric counterpart
        }
}

}

namespace matching {

namespace impl {
kronecker_approximation::kronecker_approximation(const matrix_t& A, const matrix_t& B, int max_degree) {
    // first element is 1 and it's constant
    --max_degree;
    
    m = A.cols();
    for (int t = 0; t < max_degree; t++) {
        // zero matrix of necessary size
        Matrix outer_sum = Matrix::Zero(A.cols() * A.cols());
        auto original = outer_sum;
        for (int i = 0; i < t; i++) {
            outer_sum = Eigen::KroneckerProduct(outer_sum, original);
        }

        // sum of outer products of distances
        for (int i = 0; i < A.rows(); i++) {
            for (int j = 0; j < B.rows(); j++) {
                Eigen::VectorXd d = A(i, all) - B(j, all);
                matrix_t term = d.transpose() * d;  
                auto original_term = term;
                for (int k = 0; k < t; k++) {
                    term = Eigen::KroneckerProduct(term ,original_term);
                }
                outer_sum = outer_sum + term;
            }
        }
        degrees.emplace_back(outer_sum, m, t + 1);
    }
}

// Decompose the outer product with initial dimensionality m*m at k's Kronecker degree.
// This includes logic with nkp remainders, where the nkp itself goes one level down.
// Also includes a stop condition for the approximation.
one_degree_approximation::one_degree_approximation(const matrix_t& outer, int m, int k) :m(m), k(k) {
    double init_frob = outer.squaredNorm();
    double tol = 1e-6;

    matrix_t sum_of_effects = matrix_t::Zero(outer.rows(), outer.cols());
    matrix_t residual = outer;
    while (true) {
        auto spot = one_spot_approximation(outer, m, k);
        spots.push_back(spot);

        residual = residual - spot.get_approximation();
        double curr_frob = residual.squaredNorm();
        if (curr_frob / init_frob < tol) {
            break;
        }
    }
}

one_spot_approximation::one_spot_approximation(const matrix_t& outer, int m, int k) {
    // k means kronecker's degree, m - dimensions of the matrix
    int dim = outer.cols();
    matrix_t resid = outer;
    for (int i = 0; i < k - 1; i++) {
        dim = dim / m;   
        auto result = nkp(resid, m, m, dim, dim);
        matrices.push_back(result.B);
        resid = result.C;
    }
    matrices.push_back(resid);
}

double one_spot_approximation::calculate(const matrix_t& sigma, double c_e) const {
    double res = 1.0;
    for (const auto& matrix: matrices) {
        int dim = matrix.cols();
        double trace = 0.0;
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                trace += c_e * sigma(i, j) * matrix(j, i);
            }
        }
        res *= trace;
    }
    return res;
}

one_spot_approximation::one_spot_approximation(int k, int m) :k(k), m(m) {}

one_degree_approximation::one_degree_approximation(int m, int k) :m(m), k(k) {}

kronecker_approximation::kronecker_approximation() {}

std::ostream& operator<<(std::ostream& os, const one_spot_approximation& spot) {
    if (!(os.flags() & std::ios::binary)) {
        throw std::runtime_error("Output stream must be in binary mode");
    }

    os.write("SPOT", 4);

    uint32_t n = spot.matrices.size();
    os.write(reinterpret_cast<const char*>(&n), sizeof(n));

    uint32_t k = spot.k;
    os.write(reinterpret_cast<const char*>(&k), sizeof(k));

    for (const auto& matrix: spot.matrices) {
        os.write("MATX", 4);
        write_binary(os, matrix);
    }
}

std::istream& operator>>(std::istream& is, one_spot_approximation& spot) {
    if (!(is.flags() & std::ios::binary)) {
        throw std::runtime_error("Input stream must be in binary mode");
    }

    char buf[4];
    is.read(buf, 4);
    if (std::string(buf, 4) != "SPOT") {
        throw std::runtime_error("Invalid file format");
    }

    uint32_t n;
    is.read(reinterpret_cast<char*>(&n), sizeof(n));

    spot.matrices.clear();
    for (size_t i = 0; i < n; i++) {
        is.read(buf, 4);
        if (std::string(buf, 4) != "MATX") {
            throw std::runtime_error("Invalid file format");
        }

        matrix_t matrix;
        read_binary(is, matrix);

        if (matrix.rows() != matrix.cols() || matrix.rows() != spot.m) {
            throw std::runtime_error("Invalid matrix size");
        }

        spot.matrices.push_back(matrix);
    }

    if(spot.matrices.size() != spot.k) {
        throw std::runtime_error("Invalid number of matrices");
    }
}

std::ostream& operator<<(std::ostream& os, const one_degree_approximation& spot) {
    if (!(os.flags() & std::ios::binary)) {
        throw std::runtime_error("Output stream must be in binary mode");
    }

    os.write("DEGR", 4);

    // write the number of spots
    uint32_t n = spot.spots.size();
    os.write(reinterpret_cast<const char*>(&n), sizeof(n));

    for (const auto& spot: spot.spots) {
        os << spot;
    }
}

double one_degree_approximation::calculate(const matrix_t& sigma, double c_e) const {
    double res = 0.0;
    for (const auto& spot: spots) {
        res += spot.calculate(sigma, c_e);
    }
    return res;
}

std::istream& operator>>(std::istream& is, one_degree_approximation& obj) {
    if (!(is.flags() & std::ios::binary)) {
        throw std::runtime_error("Input stream must be in binary mode");
    }

    char buf[4];
    is.read(buf, 4);
    if (std::string(buf, 4) != "DEGR") {
        throw std::runtime_error("Invalid file format");
    }

    uint32_t n;
    is.read(reinterpret_cast<char*>(&n), sizeof(n));

    obj.spots.clear();
    for (size_t i = 0; i < n; i++) {
        one_spot_approximation spot(obj.m, obj.k);
        is >> spot;
        obj.spots.push_back(spot);
    }
}

std::ostream& operator<<(std::ostream& os, const kronecker_approximation& spot) {
    if (!(os.flags() & std::ios::binary)) {
        throw std::runtime_error("Output stream must be in binary mode");
    }

    os.write("KRON", 4);

    // write a constant number to check endianess
    uint32_t magick = 0x12345678;
    os.write(reinterpret_cast<const char*>(&magick), sizeof(magick));

    // write the number of degrees
    uint32_t n = spot.degrees.size();
    os.write(reinterpret_cast<const char*>(&n), sizeof(n));

    uint32_t m = spot.m;
    os.write(reinterpret_cast<const char*>(&m), sizeof(m));

    for (const auto& degree: spot.degrees) {
        os << degree;
    }
}

double kronecker_approximation::calculate(const matrix_t& sigma, double c_e) const {
    double res = 0.0;
    double factorial = 1.0;
    for (int i = 0; i < degrees.size(); i++) {
        double deg_res =  degrees[i].calculate(sigma, c_e);
        factorial /= (i + 1);
        res += deg_res * factorial;
    }
    return res;
}

std::istream& operator>>(std::istream& is, kronecker_approximation& obj) {
    if (!(is.flags() & std::ios::binary)) {
        throw std::runtime_error("Input stream must be in binary mode");
    }

    char buf[4];
    is.read(buf, 4);
    if (std::string(buf, 4) != "KRON") {
        throw std::runtime_error("Invalid file format");
    }

    // checking endianess
    uint32_t magick;
    is.read(reinterpret_cast<char*>(&magick), sizeof(magick));
    if (magick != 0x12345678) {
        throw std::runtime_error("Invalid file format. Endianess mismatch.");
    }

    uint32_t n;
    is.read(reinterpret_cast<char*>(&n), sizeof(n));

    uint32_t m;
    is.read(reinterpret_cast<char*>(&m), sizeof(m));

    obj.degrees.clear();
    for (size_t i = 0; i < n; i++) {
        one_degree_approximation degree(m, i + 1);
        is >> degree;
        obj.degrees.push_back(degree);
    }
}

matrix_t one_spot_approximation::get_approximation() const {
    matrix_t result = matrices[0];
    for (size_t i = 1; i < matrices.size(); i++) {
        result = Eigen::KroneckerProduct(result, matrices[i]);
    }
    return result;
}

}

kronecker_calculator::kronecker_calculator(std::string filename) :fin(filename, std::ios::binary) {
    if (!fin) {
        throw std::runtime_error("Cannot open file for reading");
    }

    char buf[5];
    fin.read(buf, 5);
    if (std::string(buf, 5) != "SCORE") {
        throw std::runtime_error("Invalid file format");
    }    
}

kronecker_preprocessor::kronecker_preprocessor(std::shared_ptr<Matrix> _matrix, 
                                               std::vector<std::vector<int>> _clusters, 
                                               std::string filename) :matrix(_matrix), clusters(_clusters), fout(filename, std::ios::binary) {
    if (!fout) {
        throw std::runtime_error("Cannot open file for writing");
    }

    fout.write("SCORE", 5);

    // write the number of clusters
    uint32_t n = clusters.size();
    fout.write(reinterpret_cast<const char*>(&n), sizeof(n));
}

std::vector<std::vector<double>> kronecker_calculator::calculate(Eigen::MatrixXd sigma, double c_e) {
    uint32_t n;
    fin.read(reinterpret_cast<char*>(&n), sizeof(n));

    std::vector<std::vector<double>> result(n, std::vector<double>(n, 0)); 

    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            impl::kronecker_approximation approx;
            fin >> approx;
            double value = approx.calculate(sigma, c_e);
            if (i != j) {
                value *= 0.5;
            }
            result[i][j] = value;
            result[j][i] = value;
        }
    }

    return result;
}


void kronecker_preprocessor::process(unsigned threads, unsigned batch_size, unsigned max_degree) {
    cxxpool::thread_pool pool(threads);
    
    auto n = clusters.size();
    int curr_in_batch = 0;
    for (size_t i = 0; i < n; i++) {
        Rcpp::checkUserInterrupt();
        // Threading
        for (size_t j = i; j < n; j++) {
            futures.push_back(pool.push([this, i, j, max_degree]() -> std::unique_ptr<impl::kronecker_approximation> {
                int k = matrix->cols();
                return std::make_unique<impl::kronecker_approximation>((*matrix)(clusters[i], all), (*matrix)(clusters[j], all), max_degree);
            }));
            if (futures.size() == batch_size || (i == n - 1 && j == n - 1)) {
                write_futures();
            }
        }
    }
}

void kronecker_preprocessor::write_futures() {
    for (auto& future: futures) {
        auto approx = future.get();
        fout << *approx;
    }
    futures.clear();
}       

}