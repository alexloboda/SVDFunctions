#ifndef SRC_HW_H
#define SRC_HW_H

#include "matching.h"
#include <functional>
#include <vector>

namespace matching {

double chi2_aux(double obs, double exp);

bool check_counts(unsigned hom_ref, unsigned het, unsigned hom, double maf = 0.05, int mac = 10,
                  double chi2_bdry = 10.82757);

std::vector<bool> check_user_counts(const std::vector<Counts>& case_counts);

double chi2(int k, int n, const std::function<double(double)>& qchisq);

}

#endif //SRC_HW_H
