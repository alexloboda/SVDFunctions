#ifndef MAHALANOBIS_DISTANCES_H
#define MAHALANOBIS_DISTANCES_H

#include <vector>
#include <memory>
#include <RcppEigen.h>

namespace mvn {

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

class mahalanobis_distances {
protected:
    std::vector<double> dist;
    std::vector<std::vector<double>> inter;
    std::vector<double> diag;
    std::shared_ptr<const Matrix> X;
    Matrix S_inv;
    bool calc_interpoint;

public:
    mahalanobis_distances(std::shared_ptr<const Matrix> X, const Matrix& cov, const Vector& mean);
    void calculate_interpoint();
    Matrix get_sigma() const;
    double distance(unsigned i) const;
    double interpoint_distance(unsigned i, unsigned j) const;
};

} // namespace mvn

#endif // MAHALANOBIS_DISTANCES_H
