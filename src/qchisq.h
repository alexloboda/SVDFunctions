#ifndef SRC_QCHISQ_H
#define SRC_QCHISQ_H

#include <vector>
#include <functional>

class qchi2 {
    std::vector<double> precomputed;
    double a;
public:
    qchi2(const std::vector<double>& precomputed_values);
    double call(double p);
    std::function<double(double)> function();
};


#endif //SRC_QCHISQ_H
