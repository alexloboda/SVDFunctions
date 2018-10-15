// [[Rcpp::depends(BH)]]

#include <cmath>
#include <boost/math/special_functions/beta.hpp>

#include "lm.h"

double pval_t(double t, double df) {
    if (df > 50) {
        return 1 + boost::math::erf(-std::abs(t) / sqrt(2));
    }
    double x = df / (t * t + df);
    return boost::math::ibeta(df / 2, 0.5, x);
}

lm::lm(unsigned _n) :lt(_n * 2), n(_n), b(_n), R(4), v(_n), vM(_n), Q(_n * 2) {}

double lm::norm(const std::vector<double>& v) {
    double value = 0.0;
    for (double x: v) {
        value += x * x;
    }
    return std::sqrt(value);
}

void lm::set(int i, double x, double y, double w, bool bias) {
    double sq = std::sqrt(w);
    if (bias) {
        lt[i] = sq;
    }
    lt[i + n] = x * sq;
    b[i] = y * sq;
}

void lm::solve() {
    preprocess(v, lt, 0, n, 0);

    // Q_1 * A
    std::vector<double> RM(n);
    double cross = 0.0;
    for (int j = 0; j < n; j++) {
        cross += v[j] * lt[n + j];
    }
    for (int i = 1; i < n; i++) {
        RM[i] = lt[n + i] - 2 * v[i] * cross;
    }

    vM[0] = 0;
    preprocess(vM, RM, 0, n, 1);
    compute_Q();
    compute_R();
    solve_system();
}

void lm::preprocess(std::vector<double>& ret, const std::vector<double>& lt, int l, int r, int s) {
    int n = r - l;
    for (int i = s; i < n; i++) {
        ret[i] = lt[i + l];
    }
    ret[s] -= norm(ret);
    double l2 = norm(ret);
    if (l2 > 0) {
        for (int i = 0; i < n; i++) {
            ret[i] /= l2;
        }
    }
}

void lm::compute_Q() {
    double cross = 0.0;
    for (int k = 0; k < n; k++) {
        cross += v[k] * vM[k];
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 2; j++) {
            double sum = 0;
            sum += 4 * v[i] * vM[j] * cross;
            Q[i + j * n] = sum - 2 * (v[i] * v[j] + vM[i] * vM[j]);
        }
    }
    Q[0] += 1;
    Q[n + 1] += 1;
}

void lm::compute_R() {
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            if (i == 1 && j == 0) {
                continue;
            }
            double sum = 0;
            for (int k = 0; k < n; k++) {
                sum += Q[k + i * n] * lt[k + j * n];
            }
            R[i + j * 2] = sum;
        }
    }
    R[1] = 0;

}

void lm::solve_system() {
    double b1 = 0.0;
    double b2 = 0.0;
    for (int i = 0; i < n; i++) {
        b1 += Q[i] * b[i];
        b2 += Q[i + n] * b[i];
    }
    k = b2 / R[3];
    bias = (b1 - k * R[2]) / R[0];
}

double lm::compute_t(double df) {
    double a = R[0] * R[0];
    double b = R[0] * R[2];
    double c = R[2] * R[2] + R[3] * R[3];
    double inv =  a / (a * c - b * b);
    return k / std::sqrt((compute_rss() / df) * inv);

}

double lm::compute_rss() {
    double rss = 0;
    for (int i = 0; i < n; i++) {
        rss += pow(bias * lt[i] + k * lt[n + i] - b[i], 2.0);
    }
    return rss;
}

double lm::get_lambda() {
    return k;
}

double lm::get_bias() {
    return bias;
}

