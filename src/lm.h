#ifndef SRC_LM_H
#define SRC_LM_H

#include <vector>

double pval_t(double t, double df);

class lm {
    std::vector<double> lt;
    std::vector<double> b;
    std::vector<double> v;
    std::vector<double> vM;

    std::vector<double> R;
    std::vector<double> Q;
    double k;
    double bias;
    int n;

    double norm(const std::vector<double>& v);
    void preprocess(std::vector<double>& ret, const std::vector<double>& lt, int i, int j, int s);
    void compute_Q();
    void compute_R();
    double compute_rss();
    void solve_system();
public:
    explicit lm(unsigned n);
    void set(int i, double x, double y, double w, bool bias = true);
    void solve();
    double get_lambda();
    double get_bias();

    double compute_t(double df);
};


#endif //SRC_LM_H
