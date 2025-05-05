#include "include/mahalanobis_distances.h"
#include <stdexcept>
#include <cmath>

namespace mvn{

mahalanobis_distances::mahalanobis_distances(std::shared_ptr<const Matrix> _X, const Matrix& S, const Vector& mean)
        :dist(X->cols()), calc_interpoint(false), X(_X) {
    Eigen::FullPivHouseholderQR<Matrix> qr(S);
    if (!qr.isInvertible()) {
        throw std::logic_error("Non-invertible matrix. Must not happen.");
    }

    S_inv = qr.inverse();

    Vector ximu = X->transpose() * S_inv * mean;
    Vector muxi = mean.transpose() * S_inv * *X;
    double mumu = mean.transpose() * S_inv * mean;

    diag.resize(X->cols());

    for (auto i = 0; i < X->cols(); i++) {
        diag[i] = X->col(i).transpose() * S_inv * X->col(i);    
    }

    for (auto i = 0; i < X->cols(); i++) {
        dist[i] = diag[i] - ximu(i) - muxi(i) + mumu;
    }
}

void mahalanobis_distances::calculate_interpoint() {
    Matrix distances = X->transpose() * S_inv * *X;
    inter.resize(X->cols());

    for (auto i = 0; i < X->cols(); i++) {
        inter[i].resize(X->cols());
        for (auto j = 0; j < X->cols(); j++) {
            inter[i][j] = diag[i] - 2 *  distances(i, j) + diag[j];
        }
    }

    calc_interpoint = true;
}

Matrix mahalanobis_distances::get_sigma() const {
    return S_inv;
}

double mahalanobis_distances::distance(unsigned el) const {
    return dist[el];
}

double mahalanobis_distances::interpoint_distance(unsigned i, unsigned j) const {
    if (!calc_interpoint) {
        throw std::logic_error("Interpoint distances not calculated");
    }
    return inter[i][j];
}

} 