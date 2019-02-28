#include "include/qchisq.h"
#include <fstream>
#include <cmath>

qchi2::qchi2(const std::vector<double>& precomputed_values) {
    precomputed = precomputed_values;
    a = precomputed.size() > 10 ? 0.5 : 3.0 / 8;
}

double qchi2::call(double p) {
    p = 1 - p;
    int k = (int)std::lround(p * (precomputed.size() + 1 - 2 * a) + a) - 1;
    k = std::max(k, 0);
    k = std::min(k, (int)precomputed.size() - 1);
    return precomputed[k];
}

std::function<double(double)> qchi2::function() {
    return [this](double p) {return call(p);};
}
