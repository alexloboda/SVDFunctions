#include "include/kronecker.h"
#include <limits>
#include <cmath>
#include <random>
#include <iomanip>

// IRLBA headers
#include "include/third-party/irlba/compute.hpp"
#include "include/third-party/irlba/Options.hpp"

namespace matching {
namespace impl {

// Integer power helper for stability and speed
static inline double ipow(double x, int p) {
    double r = 1.0;
    for (int i = 0; i < p; ++i) r *= x;
    return r;
}

// one_spot_approximation
one_spot_approximation::one_spot_approximation(int m, int k) :k(k), m(m), v_(m), lambda_(0.0) {}

one_spot_approximation::one_spot_approximation(const vector_t& v, double lambda, int m, int k)
    : k(k), m(m), v_(v), lambda_(lambda) {
    double nrm = v_.norm();
    if (nrm > 0) v_ /= nrm;
    KRON_DLOG("  SPOT created: k=" << k << ", m=" << m << ", |v|=" << v_.norm() << ", lambda=" << lambda_);
}

double one_spot_approximation::calculate(const matrix_t& sigma, double c_e) const {
    double q = v_.transpose() * sigma * v_;
    double contrib = lambda_ * ipow(c_e * q, k);
    KRON_DLOG("    Spot: k=" << k << ", q=" << q << ", contrib=" << contrib);
    return contrib;
}

one_degree_approximation::one_degree_approximation(int m, int k) :m(m), k(k) {}

int one_degree_approximation::n_spots() const { return static_cast<int>(spots.size()); }

double one_degree_approximation::calculate(const matrix_t& sigma, double c_e) const {
    double res = 0.0;
    for (const auto& spot: spots) res += spot.calculate(sigma, c_e);
    return res;
}

kronecker_approximation::kronecker_approximation() {}

kronecker_approximation::kronecker_approximation(const matrix_t& A, const matrix_t& B, int max_degree)
    : kronecker_approximation(A, B, max_degree, tpm_params{}) {}

kronecker_approximation::kronecker_approximation(const matrix_t& A, const matrix_t& B, int max_degree, const tpm_params& params) {
    // Build degree-wise components using stochastic Tensor Power Method (TPM)
    m = A.cols();
    const int d = m;
    const int Dmax = max_degree;

    const int nA = A.rows();
    const int nB = B.rows();
    const std::size_t total_pairs = static_cast<std::size_t>(nA) * static_cast<std::size_t>(nB);

    const std::size_t cap = params.sample_cap;
    const std::size_t S = std::min<std::size_t>(total_pairs, cap);
    const int max_components = params.max_components;
    const int max_iters = params.max_iters;
    const double tol = params.tol;

    std::uint64_t seed = params.seed;
    if (seed == 0) {
        seed = 0x9e3779b97f4a7c15ULL;
        seed ^= static_cast<std::uint64_t>(nA) + 0x9e3779b97f4a7c15ULL + (seed<<6) + (seed>>2);
        seed ^= static_cast<std::uint64_t>(nB) + 0x9e3779b97f4a7c15ULL + (seed<<6) + (seed>>2);
        seed ^= static_cast<std::uint64_t>(d)  + 0x9e3779b97f4a7c15ULL + (seed<<6) + (seed>>2);
    }
    KRON_DLOG("TPM init: m=" << m << ", nA=" << nA << ", nB=" << nB
              << ", Dmax=" << Dmax << ", cap=" << cap << ", S=" << S
              << ", max_components=" << max_components << ", max_iters=" << max_iters
              << ", tol=" << tol << ", seed=" << seed);

    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<int> distA(0, nA ? nA - 1 : 0);
    std::uniform_int_distribution<int> distB(0, nB ? nB - 1 : 0);
    std::normal_distribution<double> gauss(0.0, 1.0);

    // Pre-sample S pairs of diffs
    std::vector<vector_t> diffs;
    diffs.reserve(S);
    if (nA > 0 && nB > 0) {
        for (std::size_t s = 0; s < S; ++s) {
            int ia = distA(rng);
            int jb = distB(rng);
            diffs.emplace_back((A.row(ia) - B.row(jb)).transpose());
        }
    }

    const double scale = (S > 0 && S < total_pairs) ? static_cast<double>(total_pairs) / static_cast<double>(S) : 1.0;
    KRON_DLOG("Sampling done: total_pairs=" << total_pairs << ", S=" << S << ", scale=" << scale);

    // Rescale diffs by a uniform factor r; compensate downstream
    double r = 1.0;
    if (!diffs.empty()) {
        double mean_sq = 0.0;
        for (const auto& dvec : diffs) mean_sq += dvec.squaredNorm();
        mean_sq /= static_cast<double>(diffs.size());
        if (mean_sq > 0.0 && std::isfinite(mean_sq)) r = 1.0 / std::sqrt(mean_sq);
        for (auto& dvec : diffs) dvec *= r;
    }
    const double inv_r  = (r != 0.0) ? (1.0 / r) : 1.0;
    const double inv_r2 = (r != 0.0) ? (1.0 / (r * r)) : 1.0;
    KRON_DLOG("Rescale diffs: r=" << r << ", inv_r=" << inv_r << ", inv_r2=" << inv_r2);

    degrees.clear();
    double inv_fact = 1.0;
    double inv_r2_pow = 1.0;
    for (int deg = 1; deg <= Dmax; ++deg) {
        inv_r2_pow *= inv_r2; // 1 / r^{2*deg}
        const double comp_deg = inv_r2_pow;
        KRON_DLOG("Degree " << deg << " start, inv_fact=" << std::scientific << inv_fact << ", comp_deg=" << comp_deg);
        one_degree_approximation degree_obj(m, deg);
        if (diffs.empty()) {
            degrees.push_back(degree_obj);
            inv_fact /= (deg + 1);
            if (inv_fact == 0.0 || !std::isfinite(inv_fact)) break;
            continue;
        }

        double true_scalar_I_deg = 0.0;
#if KRON_DEBUG
        {
            for (const auto& diff : diffs) true_scalar_I_deg += ipow(diff.squaredNorm(), deg);
            true_scalar_I_deg *= (scale * inv_fact) * comp_deg;
            KRON_DLOG("  Degree " << deg << " baseline T_k(I)=" << true_scalar_I_deg);
        }
#endif

        std::vector<vector_t> prev_w;
        std::vector<double> prev_lambda;
        prev_w.reserve(params.max_components);
        prev_lambda.reserve(params.max_components);

        // Store raw lambdas prior to ALS refit for diagnostics
        std::vector<double> raw_lambdas; raw_lambdas.reserve(params.max_components);

        double cum_lambda = 0.0;
        for (int comp = 0; comp < max_components; ++comp) {
            const int small_lambda_retry_limit = 5;
            int small_lambda_tries_left = small_lambda_retry_limit;
            bool saved_component = false;
            while (true) {
                vector_t w(m);
                for (int t = 0; t < m; ++t) w[t] = gauss(rng);
                if (w.norm() == 0) w.setConstant(1.0 / std::sqrt(static_cast<double>(m))); else w.normalize();

                int tries_left = params.restart_tries;
                while (true) {
                    vector_t w_old = w;

                    if (params.use_irls) {
                        int inner = std::max(1, params.irls_inner_iters);
                        for (int ir = 0; ir < inner; ++ir) {
                            const bool use_pre = params.irls_pre_deflate && !prev_w.empty();
                            // If experimental pre-deflation enabled, build residuals with high-order aware scaling:
                            // r_s = d_s - sum_j prev_lambda_j^{1/(2*deg)} (v_j·d_s) v_j (heuristic to reduce (v_j·r_s)^{2k})
                            std::vector<vector_t> residuals;
                            if (use_pre) {
                                residuals.reserve(S);
                                for (std::size_t s = 0; s < S; ++s) {
                                    vector_t r_s = diffs[s];
                                    for (std::size_t j = 0; j < prev_w.size(); ++j) {
                                        double proj = prev_w[j].dot(r_s);
                                        double scale_j = 0.0;
                                        if (prev_lambda[j] > 0.0) scale_j = std::pow(prev_lambda[j], 1.0 / (2.0 * deg)); // positive part
                                        r_s.noalias() -= scale_j * proj * prev_w[j];
                                    }
                                    residuals.emplace_back(std::move(r_s));
                                }
                            }
                            std::vector<double> weights(S);
                            double maxw=0.0,minw=std::numeric_limits<double>::infinity();
                            std::size_t n_floor=0,n_cap=0; double sumw=0.0,sumw2=0.0; double max_abs_proj=0.0;
                            for (std::size_t s=0;s<S;++s){
                                const vector_t &vec_s = use_pre? residuals[s]: diffs[s];
                                double sdot = vec_s.dot(w);
                                double absdot = std::abs(sdot);
                                if (absdot>max_abs_proj) max_abs_proj=absdot;
                                double ws = ipow(absdot, std::max(0,2*deg-2));
                                if (!std::isfinite(ws)) ws=0.0;
                                if (ws < params.irls_eps) { ws = params.irls_eps; ++n_floor; }
                                if (params.irls_cap>0.0 && ws>params.irls_cap) { ws=params.irls_cap; ++n_cap; }
                                weights[s]=ws; sumw+=ws; sumw2+=ws*ws; if(ws>maxw)maxw=ws; if(ws<minw)minw=ws;
                            }
                            double ess = (sumw2>0.0)? (sumw*sumw/sumw2):0.0;
                            double meanw = (S? sumw/static_cast<double>(S):0.0);
                            double varw = (S? (sumw2/static_cast<double>(S))-meanw*meanw:0.0);
                            double stdw = varw>0? std::sqrt(std::max(0.0,varw)):0.0;
                            double scale_w = (maxw>0.0 && std::isfinite(maxw))? std::sqrt(maxw):1.0;
                            Eigen::MatrixXd X(static_cast<int>(S), m);
                            for (std::size_t s=0;s<S;++s){
                                double rw = (weights[s]>0.0)? std::sqrt(weights[s])/scale_w:0.0;
                                const vector_t &vec_s = use_pre? residuals[s]: diffs[s];
                                X.row(static_cast<int>(s)) = rw * vec_s.transpose();
                            }
                            Eigen::MatrixXd U,V; Eigen::VectorXd D; irlba::Options opts; opts.max_iterations=std::max(10,max_iters); opts.convergence_tolerance=std::max(1e-7,tol); opts.cap_number=true; auto ok_iter=irlba::compute(X,1,U,V,D,opts);
#if KRON_DEBUG
                            if (use_pre){ double numE=0.0,denE=0.0; for(std::size_t s=0;s<S;++s){ numE+=residuals[s].squaredNorm(); denE+=diffs[s].squaredNorm(); } double relE = (denE>0? numE/denE:0.0); KRON_DLOG("  IRLS pre-defl energy_ratio="<<relE); }
                            KRON_DLOG("  IRLS[deg="<<deg<<", comp="<<comp<<", inner="<<(ir+1)<<"/"<<inner<<"]: IRLBA restarts="<<ok_iter.second<<", converged="<<(ok_iter.first?1:0)<<", singval="<<(D.size()>0?D[0]:-1));
                            KRON_DLOG("    IRLS weights: max|proj|="<<max_abs_proj<<", w_min="<<minw<<", w_max="<<maxw<<", mean="<<meanw<<", std="<<stdw<<", ESS="<<ess<<", floored="<<n_floor<<", capped="<<n_cap<<(use_pre?" (pre-defl)":""));
#endif
                            if (V.cols()>0){
                                w = V.col(0); double nv=w.norm(); if(nv>0 && std::isfinite(nv)) w/=nv;
                                // High-order post-SVD deflation (original model) to approximate tensor residual objective
                                if(!prev_w.empty()){
                                    vector_t w_defl = w;
                                    for(std::size_t pi=0; pi<prev_w.size(); ++pi){
                                        double dotpw = prev_w[pi].dot(w_defl);
                                        double c = ipow(dotpw, std::max(0, 2*deg-1));
                                        if(std::isfinite(c)) w_defl.noalias() -= prev_lambda[pi]*c*prev_w[pi];
                                    }
                                    double nd = w_defl.norm(); if(nd>0 && std::isfinite(nd)) w = w_defl/nd;
                                }
#if KRON_DEBUG
                                double obj_val=0.0; for(const auto &diff: diffs){ double dp=diff.dot(w); obj_val += ipow(dp, 2*deg);} obj_val *= (scale*inv_fact)*comp_deg; KRON_DLOG("    IRLS post-step objective approx="<<obj_val);
#endif
                            }
                            // Early abort IRLS inner loop if ESS too low
                            if (ess < 0.02 * static_cast<double>(S)) { KRON_DLOG("    IRLS abort: ESS too low, fallback after inner loop"); break; }
                        }
                    } else {
                        // Original gradient-like TPM step
                        int iters_done = 0;
                        for (int it = 0; it < max_iters; ++it) {
                            ++iters_done;
                            vector_t v_sample = vector_t::Zero(m);
                            for (std::size_t s = 0; s < S; ++s) {
                                double sdot = diffs[s].dot(w);
                                double s_pow_m1 = ipow(sdot, 2 * deg - 1);
                                v_sample.noalias() += s_pow_m1 * diffs[s];
                            }
                            const double alpha_obj = (scale * inv_fact) * comp_deg;
                            vector_t v = alpha_obj * v_sample;
                            for (std::size_t pi = 0; pi < prev_w.size(); ++pi) {
                                double c = ipow(prev_w[pi].dot(w), 2 * deg - 1);
                                v.noalias() -= prev_lambda[pi] * c * prev_w[pi];
                            }
                            if (params.use_gradient) v *= static_cast<double>(deg);
                            double nv = v.norm();
                            if (!(nv > 0) || !std::isfinite(nv)) break;
                            w = v / nv;
                            double delta = (w - w_old).norm();
                            if (delta < tol) break;
                            w_old = w;
                        }
#if KRON_DEBUG
                        KRON_DLOG("  TPM[deg=" << deg << ", comp=" << comp << "]: iters=" << iters_done);
#endif
                    }

                    double base = 0.0; for (std::size_t s=0; s<S; ++s) base += ipow(diffs[s].dot(w), 2*deg); base *= (scale*inv_fact)*comp_deg;
                    double defl = 0.0; for (std::size_t pi=0; pi<prev_w.size(); ++pi) defl += prev_lambda[pi]*ipow(prev_w[pi].dot(w), 2*deg);
                    double lambda = base - defl;
                    double TkI_est = 0.0; for(const auto &diff: diffs) TkI_est += ipow(diff.squaredNorm(), deg); TkI_est *= (scale*inv_fact)*comp_deg;
                    if (params.lambda_guard_factor > 0.0 && std::isfinite(lambda) && std::isfinite(TkI_est)) {
                        double guard = params.lambda_guard_factor * std::abs(TkI_est);
                        if (std::abs(lambda) > guard) {
#if KRON_DEBUG
                            KRON_DLOG("  Lambda guard triggered: |lambda|="<<std::abs(lambda)<<" > "<<guard<<" ; discarding component and switching to gradient path for this attempt");
#endif
                            if (params.use_irls) {
                                // Fallback: disable IRLS for remainder of this component attempt
                                // Re-run a gradient style single update to replace w
                                vector_t w_grad = w; for(int it=0; it< std::min(5, max_iters); ++it){ vector_t v_sample=vector_t::Zero(m); for(std::size_t s=0;s<S;++s){ double sdot=diffs[s].dot(w_grad); double s_pow_m1=ipow(sdot, 2*deg-1); v_sample.noalias() += s_pow_m1*diffs[s]; } vector_t v = (scale*inv_fact)*comp_deg * v_sample; for(std::size_t pi=0; pi<prev_w.size(); ++pi){ double c=ipow(prev_w[pi].dot(w_grad), 2*deg-1); v.noalias() -= prev_lambda[pi]*c*prev_w[pi]; } double nv=v.norm(); if(!(nv>0) || !std::isfinite(nv)) break; w_grad = v/nv; }
                                w = w_grad;
                                base = 0.0; for (std::size_t s=0; s<S; ++s) base += ipow(diffs[s].dot(w), 2*deg); base *= (scale*inv_fact)*comp_deg;
                                defl = 0.0; for (std::size_t pi=0; pi<prev_w.size(); ++pi) defl += prev_lambda[pi]*ipow(prev_w[pi].dot(w), 2*deg);
                                lambda = base - defl;
                            }
                        }
                    }
#if KRON_DEBUG
                    KRON_DLOG("  Component candidate[deg="<<deg<<", comp="<<comp<<"]: base="<<base<<", defl="<<defl<<", lambda_raw="<<lambda);
#endif

                    if (!std::isfinite(lambda)) break;
                    if (std::abs(lambda) < 1e-12) {
                        if (small_lambda_tries_left > 0) { --small_lambda_tries_left; continue; }
                        else break;
                    }

                    degree_obj.add_spot(one_spot_approximation(w, lambda, m, deg));
                    prev_w.push_back(w);
                    prev_lambda.push_back(lambda);
                    raw_lambdas.push_back(lambda);
                    cum_lambda += lambda;
                    saved_component = true;
                    break;
                }
                if (!saved_component) break;
            }

        // ALS step: refit lambdas given fixed spot vectors to minimize moment mismatch along spot directions
        if (degree_obj.n_spots() > 0) {
            // Capture pre-ALS lambda summary
#if KRON_DEBUG
            {
                double sum_raw = 0.0, sum_abs_raw = 0.0, max_abs_raw = 0.0;
                for (double L : raw_lambdas) { sum_raw += L; double a = std::abs(L); sum_abs_raw += a; if (a > max_abs_raw) max_abs_raw = a; }
                double mean_abs_raw = raw_lambdas.empty()?0.0:sum_abs_raw / raw_lambdas.size();
                KRON_DLOG("  Degree " << deg << " pre-ALS: components=" << raw_lambdas.size() << ", sum_lambda=" << sum_raw
                          << ", mean_abs_lambda=" << mean_abs_raw << ", max_abs_lambda=" << max_abs_raw);
            }
#endif
            const int p = degree_obj.n_spots();
            Eigen::MatrixXd M(p, p);
            Eigen::VectorXd b(p);
            // Build Gram-like matrix M_ij = (v_j^T v_i)^{2k}
            int row = 0;
            std::vector<vector_t> V; V.reserve(p);
            for (const auto& si : degree_obj) V.push_back(si.v());
            for (int i = 0; i < p; ++i) {
                for (int j = 0; j < p; ++j) {
                    double cij = V[j].dot(V[i]);
                    M(i, j) = ipow(cij, 2 * deg);
                }
            }
            // Right-hand side b_i = (scale*inv_fact)*comp_deg * sum_s (d_s·v_i)^{2k}
            for (int i = 0; i < p; ++i) {
                double base_i = 0.0;
                for (const auto& diff : diffs) base_i += ipow(diff.dot(V[i]), 2 * deg);
                b[i] = (scale * inv_fact) * comp_deg * base_i;
            }
            // Tiny ridge for stability (does not bias identity-scalar much)
            const double ridge = 1e-12;
            M.diagonal().array() += ridge;
            Eigen::VectorXd lambda_new;
            bool solved = false;
            // Prefer LDLT if positive definite-ish; else fall back to QR
            Eigen::LDLT<Eigen::MatrixXd> ldlt(M);
            if (ldlt.info() == Eigen::Success) {
                lambda_new = ldlt.solve(b);
                if (ldlt.info() == Eigen::Success && lambda_new.allFinite()) solved = true;
            }
            if (!solved) {
                lambda_new = M.colPivHouseholderQr().solve(b);
                if (lambda_new.allFinite()) solved = true;
            }
            if (solved) {
                // Rebuild degree object with same v and new lambdas
                one_degree_approximation refit_obj(m, deg);
                std::vector<double> refit_lambdas; refit_lambdas.reserve(p);
                for (int i = 0; i < p; ++i) {
                    refit_obj.add_spot(one_spot_approximation(V[i], lambda_new[i], m, deg));
                    refit_lambdas.push_back(lambda_new[i]);
                }
                degree_obj = refit_obj;
#if KRON_DEBUG
                double sum_refit = 0.0, sum_abs_refit = 0.0, max_abs_refit = 0.0, max_abs_diff = 0.0, mean_abs_diff = 0.0;
                for (int i = 0; i < p; ++i) {
                    double Lr = raw_lambdas[i];
                    double Ln = refit_lambdas[i];
                    sum_refit += Ln; double a = std::abs(Ln); sum_abs_refit += a; if (a > max_abs_refit) max_abs_refit = a;
                    double diff = std::abs(Lr - Ln); mean_abs_diff += diff; if (diff > max_abs_diff) max_abs_diff = diff;
                }
                if (p > 0) mean_abs_diff /= p;
                double mean_abs_refit = p? (sum_abs_refit / p) : 0.0;
                KRON_DLOG("  Degree " << deg << " ALS: refit lambdas p=" << p << ", sum_refit=" << sum_refit
                          << ", mean_abs_refit=" << mean_abs_refit << ", max_abs_refit=" << max_abs_refit
                          << ", mean_abs_diff=" << mean_abs_diff << ", max_abs_diff=" << max_abs_diff);
#endif
                KRON_DLOG("  Degree " << deg << " ALS: refit lambdas with p=" << p);
            } else {
                KRON_DLOG("  Degree " << deg << " ALS: solve failed; keeping TPM lambdas");
            }
        }

        // Optional per-degree identity-scalar normalization to correct global scaling bias
        if (degree_obj.n_spots() > 0 && params.enforce_identity_scalar) {
            double approx_scalar_I = 0.0;
            for (const auto& spot : degree_obj) approx_scalar_I += spot.lambda();
            // true_scalar_I_deg already computed with current comp_deg
            double true_scalar_I = 0.0;
            for (const auto& diff : diffs) true_scalar_I += ipow(diff.squaredNorm(), deg);
            true_scalar_I *= (scale * inv_fact) * comp_deg;
            if (std::isfinite(true_scalar_I) && std::isfinite(approx_scalar_I) && std::abs(approx_scalar_I) > 0.0) {
                double s = true_scalar_I / approx_scalar_I;
                // Rescale all lambdas uniformly
                one_degree_approximation scaled_obj(m, deg);
                for (const auto& spot : degree_obj) {
                    scaled_obj.add_spot(one_spot_approximation(spot.v(), spot.lambda() * s, m, deg));
                }
#if KRON_DEBUG
                KRON_DLOG("  Degree " << deg << " identity-scalar enforced: scale=" << s << ", sum_before=" << approx_scalar_I << ", sum_target=" << true_scalar_I);
#endif
                degree_obj = scaled_obj;
            }
        }

        degrees.push_back(degree_obj);

#if KRON_DEBUG
        // Per-degree summary aggregates prior to detailed diagnostics
        {
            const auto& deg_spots = degrees.back();
            int p = deg_spots.n_spots();
            double sum_lambda = 0.0, sum_abs_lambda = 0.0, max_abs_lambda = 0.0;
            for (const auto& spot : deg_spots) { double L = spot.lambda(); sum_lambda += L; double a = std::abs(L); sum_abs_lambda += a; if (a > max_abs_lambda) max_abs_lambda = a; }
            double mean_abs_lambda = p? (sum_abs_lambda / p) : 0.0;
            double rel_err_I = 0.0;
            if (true_scalar_I_deg != 0.0) rel_err_I = std::abs(sum_lambda - true_scalar_I_deg) / std::abs(true_scalar_I_deg);
            KRON_DLOG("  Degree " << deg << " summary: components=" << p << ", sum_lambda=" << sum_lambda
                      << ", T_k(I)=" << true_scalar_I_deg << ", rel_err_I=" << rel_err_I
                      << ", mean_abs_lambda=" << mean_abs_lambda << ", max_abs_lambda=" << max_abs_lambda);
        }
        KRON_DLOG("Degree " << deg << " diagnostics start");
        if (!diffs.empty()) {
            if (deg == 1) {
                Eigen::MatrixXd true_M = Eigen::MatrixXd::Zero(m, m);
                for (const auto& diff : diffs) true_M.noalias() += diff * diff.transpose();
                true_M *= (scale * inv_fact) * comp_deg;

                Eigen::MatrixXd approx_M = Eigen::MatrixXd::Zero(m, m);
                for (const auto& spot : degrees.back()) approx_M.noalias() += spot.lambda() * (spot.v() * spot.v().transpose());
                double frob_M = (true_M - approx_M).norm();
                double tr_true = true_M.trace();
                double tr_approx = approx_M.trace();
                double tr_abs_err = std::abs(tr_true - tr_approx);
                double tr_rel_err = (std::abs(tr_true) > 0.0) ? (tr_abs_err / std::abs(tr_true)) : tr_abs_err;
                KRON_DLOG("  Degree 1 matrix diagnostic: Frobenius=" << frob_M << ", trace_true=" << tr_true << ", trace_approx=" << tr_approx << ", trace_rel_err=" << tr_rel_err);

                vector_t true_vec = vector_t::Zero(m);
                for (const auto& diff : diffs) true_vec.noalias() += diff;
                true_vec *= (scale * inv_fact) * inv_r;

                vector_t approx_vec = vector_t::Zero(m);
                for (const auto& spot : degrees.back()) approx_vec.noalias() += spot.lambda() * spot.v();
                const double true_n = true_vec.norm();
                const double approx_n = approx_vec.norm();
                const double l2 = (true_vec - approx_vec).norm();
                double cosang = 0.0;
                if (true_n > 0 && approx_n > 0) cosang = true_vec.dot(approx_vec) / (true_n * approx_n);
                const double rel = (true_n > 0) ? (l2 / true_n) : l2;
                KRON_DLOG("  Degree 1 vector diagnostic: L2=" << l2 << ", ||true||=" << true_n << ", ||approx||=" << approx_n << ", cos(theta)=" << cosang << ", rel_err=" << rel);
            } else if (deg == 2) {
                Eigen::MatrixXd true_tensor = Eigen::MatrixXd::Zero(m, m);
                for (const auto& diff : diffs) true_tensor.noalias() += diff * diff.transpose();
                // Use comp_deg (r^{-2k}) consistently for degree-2 instead of inv_r2
                true_tensor *= (scale * inv_fact) * comp_deg;

                Eigen::MatrixXd approx_tensor = Eigen::MatrixXd::Zero(m, m);
                for (const auto& spot : degrees.back()) approx_tensor.noalias() += spot.lambda() * (spot.v() * spot.v().transpose());
                double frob_norm_diff = (true_tensor - approx_tensor).norm();
                KRON_DLOG("  Degree 2 Frobenius norm of difference: " << frob_norm_diff);
            } else {
                double max_abs_err = 0.0;
                double mean_abs_err = 0.0;
                int cnt = 0;
                const auto& deg_spots = degrees.back();
                for (const auto& si : deg_spots) {
                    const vector_t& vi = si.v();
                    double base_i = 0.0;
                    for (const auto& diff : diffs) base_i += ipow(diff.dot(vi), 2 * deg);
                    base_i *= (scale * inv_fact) * comp_deg;
                    double defl_i = 0.0;
                    for (const auto& sj : deg_spots) { if (&sj == &si) continue; defl_i += sj.lambda() * ipow(sj.v().dot(vi), 2 * deg); }
                    double lambda_check = base_i - defl_i;
                    double err = std::abs(lambda_check - si.lambda());
                    max_abs_err = std::max(max_abs_err, err);
                    mean_abs_err += err;
                    ++cnt;
                }
                if (cnt > 0) mean_abs_err /= static_cast<double>(cnt);
                KRON_DLOG("  Degree " << deg << " lambda consistency: mean_abs_err=" << mean_abs_err << ", max_abs_err=" << max_abs_err);
            }

            double true_scalar_I = true_scalar_I_deg;
            double approx_scalar_I = 0.0;
            for (const auto& spot : degrees.back()) approx_scalar_I += spot.lambda();
            double abs_err_I = std::abs(true_scalar_I - approx_scalar_I);
            double denom_I = std::abs(true_scalar_I);
            double rel_err_I = (denom_I > 0.0) ? (abs_err_I / denom_I) : abs_err_I;
            KRON_DLOG("  Degree " << deg << " identity-scalar: true=" << true_scalar_I << ", approx=" << approx_scalar_I << ", abs_err=" << abs_err_I << ", rel_err=" << rel_err_I);

            std::mt19937_64 dbg_rng(0xC0FFEEULL + static_cast<std::uint64_t>(deg));
            std::normal_distribution<double> g(0.0, 1.0);
            const int U = 3;
            for (int ui = 0; ui < U; ++ui) {
                vector_t u(m);
                for (int t = 0; t < m; ++t) u[t] = g(dbg_rng);
                double nu = u.norm();
                if (nu == 0.0 || !std::isfinite(nu)) { u.setZero(); u[0] = 1.0; nu = 1.0; }
                u /= nu;
                double true_u = 0.0;
                for (const auto& diff : diffs) true_u += ipow(diff.dot(u), 2 * deg);
                true_u *= (scale * inv_fact) * comp_deg;
                double approx_u = 0.0;
                for (const auto& spot : degrees.back()) approx_u += spot.lambda() * ipow(spot.v().dot(u), 2 * deg);
                double abs_err_u = std::abs(true_u - approx_u);
                double rel_err_u = (std::abs(true_u) > 0.0) ? (abs_err_u / std::abs(true_u)) : abs_err_u;
                KRON_DLOG("  Degree " << deg << " proj u#" << ui << ": true=" << true_u << ", approx=" << approx_u << ", abs_err=" << abs_err_u << ", rel_err=" << rel_err_u);
            }
        }
        KRON_DLOG("Degree " << deg << " diagnostics end");
#endif

        inv_fact /= (deg + 1);
        if (inv_fact == 0.0 || !std::isfinite(inv_fact)) break;
    }
}

// factorial scaling already included in lambdas

double kronecker_approximation::calculate(const matrix_t& sigma, double c_e) const {
    if (sigma.rows() != m || sigma.cols() != m) throw std::invalid_argument("sigma must be m x m");
    if (!std::isfinite(c_e)) throw std::invalid_argument("c_e must be finite");
    if (!sigma.allFinite()) throw std::invalid_argument("sigma must have all finite entries");
    KRON_DLOG("Calculate: m=" << m << ", c_e=" << c_e << ", degrees=" << degrees.size());
    double res = 0.0;
    for (const auto& degree : degrees) {
        double part = degree.calculate(sigma, c_e);
        KRON_DLOG("  Degree k=" << degree.degree() << ": contrib=" << part);
        res += part;
    }
    KRON_DLOG("Calculate: result=" << res);
    return res;
}

// Stream operators for serialization
std::ostream& operator<<(std::ostream& os, const one_spot_approximation& spot) {
    os.write("SPOT", 4);
    uint32_t k = static_cast<uint32_t>(spot.degree());
    os.write(reinterpret_cast<const char*>(&k), sizeof(k));
    uint32_t m = static_cast<uint32_t>(spot.dim());
    os.write(reinterpret_cast<const char*>(&m), sizeof(m));
    double lambda = spot.lambda();
    os.write(reinterpret_cast<const char*>(&lambda), sizeof(lambda));
    uint32_t vdim = static_cast<uint32_t>(spot.v().size());
    if (vdim != m) throw std::runtime_error("SPOT write: vector dimension mismatch with m");
    os.write(reinterpret_cast<const char*>(&vdim), sizeof(vdim));
    for (uint32_t i = 0; i < vdim; ++i) {
        double val = spot.v()[static_cast<int>(i)];
        os.write(reinterpret_cast<const char*>(&val), sizeof(val));
    }
    if (!os) throw std::runtime_error("Failed to write SPOT block");
    return os;
}

std::istream& operator>>(std::istream& is, one_spot_approximation& spot) {
    char buf[4];
    is.read(buf, 4);
    if (!is || std::string(buf, 4) != "SPOT") throw std::runtime_error("Invalid file format: missing SPOT header");
    uint32_t k; is.read(reinterpret_cast<char*>(&k), sizeof(k));
    uint32_t m; is.read(reinterpret_cast<char*>(&m), sizeof(m));
    if (!is) throw std::runtime_error("Failed to read SPOT meta");
    spot.k = static_cast<int>(k);
    spot.m = static_cast<int>(m);
    is.read(reinterpret_cast<char*>(&spot.lambda_), sizeof(spot.lambda_));
    if (!is || !std::isfinite(spot.lambda_)) throw std::runtime_error("Invalid or non-finite lambda in SPOT");
    uint32_t vdim = 0; is.read(reinterpret_cast<char*>(&vdim), sizeof(vdim));
    if (!is) throw std::runtime_error("Failed to read SPOT vector dimension");
    if (vdim != m) throw std::runtime_error("SPOT read: vector dimension mismatch with m");
    spot.v_.resize(static_cast<int>(vdim));
    for (uint32_t i = 0; i < vdim; ++i) {
        double val = 0.0; is.read(reinterpret_cast<char*>(&val), sizeof(val));
        if (!is || !std::isfinite(val)) throw std::runtime_error("Invalid or non-finite value in SPOT vector");
        spot.v_[static_cast<int>(i)] = val;
    }
    double nrm = spot.v_.norm();
    if (!(nrm > 0) || !std::isfinite(nrm)) throw std::runtime_error("Invalid SPOT vector: non-positive or non-finite norm");
    spot.v_ /= nrm;
    return is;
}

std::ostream& operator<<(std::ostream& os, const one_degree_approximation& spot) {
    os.write("DEGR", 4);
    uint32_t n = static_cast<uint32_t>(spot.spots.size());
    os.write(reinterpret_cast<const char*>(&n), sizeof(n));
    for (const auto& s: spot.spots) os << s;
    if (!os) throw std::runtime_error("Failed to write DEGR block");
    return os;
}

std::istream& operator>>(std::istream& is, one_degree_approximation& obj) {
    char buf[4];
    is.read(buf, 4);
    if (!is || std::string(buf, 4) != "DEGR") throw std::runtime_error("Invalid file format: missing DEGR header");
    uint32_t n; is.read(reinterpret_cast<char*>(&n), sizeof(n));
    if (!is) throw std::runtime_error("Failed to read number of spots in DEGR");
    obj.spots.clear(); obj.spots.reserve(n);
    for (size_t i = 0; i < n; i++) {
        one_spot_approximation spot(obj.dim(), obj.degree());
        is >> spot;
        if (spot.dim() != obj.dim() || spot.degree() != obj.degree()) throw std::runtime_error("DEGR content mismatch: SPOT m/k do not match degree");
        obj.spots.push_back(spot);
    }
    return is;
}

std::ostream& operator<<(std::ostream& os, const kronecker_approximation& spot) {
    os.write("KRON", 4);
    uint32_t magick = 0x12345679; os.write(reinterpret_cast<const char*>(&magick), sizeof(magick));
    uint32_t version = 1; os.write(reinterpret_cast<const char*>(&version), sizeof(version));
    uint32_t n = static_cast<uint32_t>(spot.degrees.size()); os.write(reinterpret_cast<const char*>(&n), sizeof(n));
    uint32_t m = static_cast<uint32_t>(spot.m); os.write(reinterpret_cast<const char*>(&m), sizeof(m));
    for (const auto& degree: spot.degrees) os << degree;
    if (!os) throw std::runtime_error("Failed to write KRON block");
    return os;
}

std::istream& operator>>(std::istream& is, kronecker_approximation& obj) {
    char buf[4]; is.read(buf, 4);
    if (!is || std::string(buf, 4) != "KRON") throw std::runtime_error("Invalid file format");
    uint32_t magick = 0; is.read(reinterpret_cast<char*>(&magick), sizeof(magick));
    if (!is) throw std::runtime_error("Failed to read KRON magic");
    uint32_t version = 0;
    if (magick == 0x12345678) {
        version = 0;
    } else if (magick == 0x12345679) {
        is.read(reinterpret_cast<char*>(&version), sizeof(version));
        if (!is) throw std::runtime_error("Failed to read KRON version");
        if (version != 1) throw std::runtime_error("Unsupported KRON version");
    } else {
        throw std::runtime_error("Invalid file format. Endianess/magic mismatch.");
    }
    uint32_t n = 0; is.read(reinterpret_cast<char*>(&n), sizeof(n));
    uint32_t m = 0; is.read(reinterpret_cast<char*>(&m), sizeof(m));
    if (!is) throw std::runtime_error("Failed to read KRON meta");
    obj.m = static_cast<int>(m);
    obj.degrees.clear(); obj.degrees.reserve(n);
    for (size_t i = 0; i < n; i++) {
        one_degree_approximation degree(obj.m, static_cast<int>(i) + 1);
        is >> degree;
        obj.degrees.push_back(degree);
    }
    return is;
}

double kronecker_approximation::compression() const {
    int spots_n = degrees.empty() ? 0 : degrees.back().n_spots();
    int degs = static_cast<int>(degrees.size());
    double expected_spots = std::pow(m * m, std::max(0, degs - 1));
    return expected_spots > 0 ? (spots_n / expected_spots) : 0.0;
}

} // namespace impl
} // namespace matching
