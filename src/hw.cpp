#include "include/hw.h"

double matching::chi2_aux(double obs, double exp) {
    double chi = std::abs(obs - exp) - 0.5;
    return (chi * chi) / exp;
}

bool
matching::check_counts(unsigned int hom_ref, unsigned int het, unsigned int hom, double maf, int mac, double chi2_bdry) {
    if (hom < hom_ref) {
        std::swap(hom, hom_ref);
    }
    double AC = 2 * hom_ref + het;
    unsigned n = hom_ref + het + hom;
    double p = AC / (2 * n);
    double q = 1 - p;
    double chi2 = chi2_aux(hom_ref, n * p * p) + chi2_aux(het, 2 * n * p * q) + chi2_aux(hom, n * q * q);
    return p > maf && AC > mac && chi2 < chi2_bdry;
}

std::vector<bool> matching::check_user_counts(const std::vector <Counts>& case_counts) {
    std::vector<bool> mask;
    for (auto& case_count : case_counts) {
        mask.push_back(check_counts(case_count[0], case_count[1], case_count[2]));
    }
    return mask;
}

double matching::chi2(int k, int n, const std::function<double(double)>& qchisq) {
    double a = n <= 10 ? 3.0 / 8.0 : 0.5;
    return qchisq((k + 1 - a) / (n + 1 - 2 * a));
}
