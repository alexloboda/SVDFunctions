#include <cmath>
#include <limits>
#include <algorithm>
#include <numeric>
#include <random>
#include <future>

#include "include/vcf_predicting_handler.h"
#include "include/genotype_predictor.h"
#include "include/vcf_parser.h"

namespace {
    using std::size_t;

    struct SymPosDefSolver {
        // Simple Cholesky decomposition for symmetric positive definite matrix.
        // Stores lower-triangular L such that A = L * L^T.
        std::vector<double> L;
        size_t n;

        explicit SymPosDefSolver(size_t n) : L(n * n, 0.0), n(n) {}

        static inline double& at(std::vector<double>& a, size_t n, size_t i, size_t j) { return a[i * n + j]; }
        static inline const double& at(const std::vector<double>& a, size_t n, size_t i, size_t j) { return a[i * n + j]; }

        bool factorize(const std::vector<double>& A) {
            // Copy A into L then factor in-place.
            L = A;
            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j <= i; j++) {
                    double sum = at(L, n, i, j);
                    for (size_t k = 0; k < j; k++) {
                        sum -= at(L, n, i, k) * at(L, n, j, k);
                    }
                    if (i == j) {
                        if (sum <= 0.0) {
                            return false;
                        }
                        at(L, n, i, j) = std::sqrt(sum);
                    } else {
                        at(L, n, i, j) = sum / at(L, n, j, j);
                    }
                }
                for (size_t j = i + 1; j < n; j++) {
                    at(L, n, i, j) = 0.0;
                }
            }
            return true;
        }

        // Solves A x = b using stored Cholesky factors.
        std::vector<double> solve(const std::vector<double>& b) const {
            std::vector<double> y(n, 0.0);
            for (size_t i = 0; i < n; i++) {
                double sum = b[i];
                for (size_t k = 0; k < i; k++) {
                    sum -= at(L, n, i, k) * y[k];
                }
                y[i] = sum / at(L, n, i, i);
            }
            std::vector<double> x(n, 0.0);
            for (size_t i = n; i-- > 0;) {
                double sum = y[i];
                for (size_t k = i + 1; k < n; k++) {
                    sum -= at(L, n, k, i) * x[k];
                }
                x[i] = sum / at(L, n, i, i);
            }
            return x;
        }
    };

    static inline int clamp_round_to_genotype(double x) {
        int r = (int)std::llround(x);
        if (r < 0) {
            return 0;
        }
        if (r > 2) {
            return 2;
        }
        return r;
    }
}

namespace vcf {
    PredictingHandler::PredictingHandler(const std::vector<std::string>& samples, GenotypeMatrixHandler& gh,
                                         int window_size_kb, int window_size, unsigned int seed)
                                         :VariantsHandler(samples), curr_chr(-1), iterator{gh},
                                          window(window_size),
                                          thread_pool(std::thread::hardware_concurrency()),
                                          random_seed(seed) {
        auto variants = gh.desired_variants();
        int halfws = window_size_kb / 2;
        for (const Variant& v : variants) {
            Chromosome chr = v.position().chromosome();
            int pos = v.position().position();
            Range haplotype{chr, pos - halfws, pos + halfws};
            ranges.insert(haplotype);
        }
    }

    bool PredictingHandler::isOfInterest(const Variant& variant) {
        Position pos = variant.position();
        if (ranges.empty()) {
            return true;
        }
        return ranges.includes(pos);
    }

    void PredictingHandler::processVariant(const Variant& variant, std::shared_ptr<AlleleVector>& alleles) {
        if (!isOfInterest(variant)) {
            return;
        }
        Position pos = variant.position();
        if (curr_chr.num() == -1) {
            curr_chr = pos.chromosome();
        }
        if (pos.chromosome() != curr_chr) {
            cleanup();
            curr_chr = pos.chromosome();
        }

        window.add(alleles, variant);
        if (window.is_full()) {
            while (iterator.dereferencable() && (*iterator).position().position() <= window.middle_point()) {
                auto dataset = window.dataset(*iterator);
                if (dataset.second.empty()) {
                    return;
                }
                Variant v = *iterator;
                fix_labels(v, dataset);
                ++iterator;
            }
        }
    }

    void PredictingHandler::cleanup() {
        for(; iterator.dereferencable(); ++iterator) {
            Variant var = *iterator;
            auto dataset = window.dataset(var);
            if (dataset.second.empty()) {
                break;
            }
            fix_labels(var, dataset);
        }
        window.clear();
    }

    void PredictingHandler::fix_labels(const Variant& variant, const std::pair<Features, Labels>& dataset) {
        bool missing = false;
        for (auto l: dataset.second) {
            if (l == MISSING) {
                missing = true;
            }
        }
        if (!missing) {
            return;
        }
        TreeBuilder tree_builder = make_tree_builder(dataset);
        RandomForest forest{tree_builder, thread_pool, /*ntrees=*/100, random_seed};

        // First pass: produce the returned genotype row (impute only missing).
        std::size_t n_observed = 0;
        std::size_t n_missing = 0;
        std::vector<size_t> observed_idx;
        observed_idx.reserve(dataset.second.size());

        std::vector<float> labels;
        labels.reserve(dataset.second.size());
        for (size_t i = 0; i < dataset.second.size(); i++) {
            AlleleType curr = dataset.second[i];
            if (curr == MISSING) {
                ++n_missing;
                std::vector<AlleleType> features;
                features.reserve(dataset.first.size());
                for (size_t j = 0; j < dataset.first.size(); j++) {
                    features.push_back(dataset.first[j][i]);
                }
                labels.push_back((float)forest.predict(features));
            } else {
                ++n_observed;
                observed_idx.push_back(i);
                labels.push_back((float)to_int(curr));
            }
        }
        iterator.set(labels);

        // RF OOB evaluation on observed genotypes (parallelized in chunks when beneficial).
        struct EvalAccum {
            double sum_abs = 0.0;
            double sum_sq = 0.0;
            std::size_t correct = 0;
            std::size_t count = 0;
        };

        auto eval_rf_range = [&](size_t begin, size_t end) -> EvalAccum {
            EvalAccum acc;
            for (size_t k = begin; k < end; k++) {
                size_t i = observed_idx[k];
                std::vector<AlleleType> features;
                features.reserve(dataset.first.size());
                for (size_t j = 0; j < dataset.first.size(); j++) {
                    features.push_back(dataset.first[j][i]);
                }
                double pred = forest.predict_oob(features, i);
                double truth = (double)to_int(dataset.second[i]);
                double err = pred - truth;
                acc.sum_abs += std::abs(err);
                acc.sum_sq += err * err;
                int rounded = clamp_round_to_genotype(pred);
                if ((double)rounded == truth) {
                    acc.correct++;
                }
                acc.count++;
            }
            return acc;
        };

        std::size_t pool_threads = 0;
        try {
            pool_threads = thread_pool.n_threads();
        } catch (...) {
            pool_threads = 0;
        }

        EvalAccum rf_acc;
        if (observed_idx.size() >= 512 && pool_threads >= 2) {
            const size_t target_tasks = std::min<std::size_t>(pool_threads * 4, 32);
            const size_t chunk = std::max<std::size_t>((observed_idx.size() + target_tasks - 1) / target_tasks, 128);
            std::vector<std::future<EvalAccum>> futures;
            for (size_t begin = 0; begin < observed_idx.size(); begin += chunk) {
                size_t end = std::min(begin + chunk, observed_idx.size());
                futures.push_back(thread_pool.push([&eval_rf_range, begin, end]() { return eval_rf_range(begin, end); }));
            }
            for (auto& f : futures) {
                auto a = f.get();
                rf_acc.sum_abs += a.sum_abs;
                rf_acc.sum_sq += a.sum_sq;
                rf_acc.correct += a.correct;
                rf_acc.count += a.count;
            }
        } else {
            rf_acc = eval_rf_range(0, observed_idx.size());
        }

        // Additional model comparisons (evaluation only; does not affect imputed labels).
        // Compute on a bounded subsample of observed samples to keep runtime manageable.
        const size_t p_total = dataset.first.size();
        const size_t p_use = std::min((size_t)30, p_total);
        const size_t train_max = 2000;
        const size_t eval_max = 1000;

        std::mt19937 rng((unsigned)(random_seed + (unsigned)variant.position().position()));
        std::shuffle(observed_idx.begin(), observed_idx.end(), rng);
        const size_t train_n = std::min(train_max, observed_idx.size());
        const size_t eval_n = std::min(eval_max, observed_idx.size());
        std::vector<size_t> train_idx(observed_idx.begin(), observed_idx.begin() + train_n);
        std::vector<size_t> eval_idx(observed_idx.begin(), observed_idx.begin() + eval_n);

        // Ridge regression LOOCV
        double ridge_mae = std::numeric_limits<double>::quiet_NaN();
        double ridge_rmse = std::numeric_limits<double>::quiet_NaN();
        double ridge_acc = std::numeric_limits<double>::quiet_NaN();
        if (p_use > 0 && train_n >= 10 && eval_n >= 10) {
            // feature means for missing imputation
            std::vector<double> feat_mean(p_use, 0.0);
            std::vector<size_t> feat_cnt(p_use, 0);
            for (size_t ti = 0; ti < train_n; ti++) {
                size_t idx = train_idx[ti];
                for (size_t j = 0; j < p_use; j++) {
                    AlleleType a = dataset.first[j][idx];
                    if (a != MISSING) {
                        feat_mean[j] += (double)to_int(a);
                        feat_cnt[j] += 1;
                    }
                }
            }
            for (size_t j = 0; j < p_use; j++) {
                feat_mean[j] = feat_cnt[j] == 0 ? 0.0 : feat_mean[j] / (double)feat_cnt[j];
            }

            // Design: intercept + p_use features
            const size_t d = p_use + 1;
            const double lambda = 1.0;
            std::vector<double> XtX(d * d, 0.0);
            std::vector<double> Xty(d, 0.0);

            auto add_outer = [&](const std::vector<double>& x, double y) {
                for (size_t a = 0; a < d; a++) {
                    Xty[a] += x[a] * y;
                    for (size_t b = 0; b < d; b++) {
                        XtX[a * d + b] += x[a] * x[b];
                    }
                }
            };

            std::vector<std::vector<double>> X_train(train_n, std::vector<double>(d, 0.0));
            std::vector<double> y_train(train_n, 0.0);
            for (size_t ti = 0; ti < train_n; ti++) {
                size_t idx = train_idx[ti];
                X_train[ti][0] = 1.0;
                for (size_t j = 0; j < p_use; j++) {
                    AlleleType a = dataset.first[j][idx];
                    X_train[ti][j + 1] = (a == MISSING) ? feat_mean[j] : (double)to_int(a);
                }
                y_train[ti] = (double)to_int(dataset.second[idx]);
                add_outer(X_train[ti], y_train[ti]);
            }

            // Ridge penalty (do not penalize intercept)
            for (size_t k = 1; k < d; k++) {
                XtX[k * d + k] += lambda;
            }

            SymPosDefSolver solver(d);
            if (solver.factorize(XtX)) {
                std::vector<double> beta = solver.solve(Xty);

                std::vector<unsigned char> in_train(dataset.second.size(), 0);
                for (size_t t = 0; t < train_n; t++) {
                    in_train[train_idx[t]] = 1;
                }

                double sum_abs_r = 0.0;
                double sum_sq_r = 0.0;
                size_t correct_r = 0;
                size_t used_r = 0;

                auto eval_ridge_range = [&](size_t begin, size_t end) -> EvalAccum {
                    EvalAccum acc;
                    for (size_t ei = begin; ei < end; ei++) {
                        size_t idx = eval_idx[ei];

                        std::vector<double> x(d, 0.0);
                        x[0] = 1.0;
                        for (size_t j = 0; j < p_use; j++) {
                            AlleleType a = dataset.first[j][idx];
                            x[j + 1] = (a == MISSING) ? feat_mean[j] : (double)to_int(a);
                        }
                        double truth = (double)to_int(dataset.second[idx]);
                        double yhat = 0.0;
                        for (size_t k = 0; k < d; k++) {
                            yhat += x[k] * beta[k];
                        }

                        double yhat_loo = yhat;
                        if (in_train[idx] == 1) {
                            std::vector<double> v = solver.solve(x);
                            double h = 0.0;
                            for (size_t k = 0; k < d; k++) {
                                h += x[k] * v[k];
                            }
                            double denom = 1.0 - h;
                            if (std::abs(denom) > 1e-8) {
                                double e = truth - yhat;
                                yhat_loo = truth - e / denom;
                            }
                        }

                        double err = yhat_loo - truth;
                        acc.sum_abs += std::abs(err);
                        acc.sum_sq += err * err;
                        int rounded = clamp_round_to_genotype(yhat_loo);
                        if ((double)rounded == truth) {
                            acc.correct++;
                        }
                        acc.count++;
                    }
                    return acc;
                };

                EvalAccum ridge_accum;
                if (eval_n >= 256 && pool_threads >= 2) {
                    const size_t target_tasks = std::min<std::size_t>(pool_threads * 4, 32);
                    const size_t chunk = std::max<std::size_t>((eval_n + target_tasks - 1) / target_tasks, 64);
                    std::vector<std::future<EvalAccum>> futures;
                    for (size_t begin = 0; begin < eval_n; begin += chunk) {
                        size_t end = std::min(begin + chunk, eval_n);
                        futures.push_back(thread_pool.push([&eval_ridge_range, begin, end]() { return eval_ridge_range(begin, end); }));
                    }
                    for (auto& f : futures) {
                        auto a = f.get();
                        ridge_accum.sum_abs += a.sum_abs;
                        ridge_accum.sum_sq += a.sum_sq;
                        ridge_accum.correct += a.correct;
                        ridge_accum.count += a.count;
                    }
                } else {
                    ridge_accum = eval_ridge_range(0, eval_n);
                }

                used_r = ridge_accum.count;
                sum_abs_r = ridge_accum.sum_abs;
                sum_sq_r = ridge_accum.sum_sq;
                correct_r = ridge_accum.correct;

                if (used_r > 0) {
                    ridge_mae = sum_abs_r / (double)used_r;
                    ridge_rmse = std::sqrt(sum_sq_r / (double)used_r);
                    ridge_acc = (double)correct_r / (double)used_r;
                }
            }
        }

        // kNN leave-one-out (within train subset); distance in genotype space over p_use features.
        double knn_mae = std::numeric_limits<double>::quiet_NaN();
        double knn_rmse = std::numeric_limits<double>::quiet_NaN();
        double knn_acc = std::numeric_limits<double>::quiet_NaN();
        if (p_use > 0 && train_n >= 20 && eval_n >= 10) {
            const size_t K = 20;
            const double eps_w = 1e-6;
            auto eval_knn_range = [&](size_t begin, size_t end) -> EvalAccum {
                EvalAccum acc;
                for (size_t ei = begin; ei < end; ei++) {
                    size_t idx = eval_idx[ei];
                    double truth = (double)to_int(dataset.second[idx]);

                    std::vector<std::pair<int, double>> best; // (distance, label)
                    best.reserve(K + 1);

                    for (size_t tj = 0; tj < train_n; tj++) {
                        size_t jdx = train_idx[tj];
                        if (jdx == idx) {
                            continue;
                        }
                        int dist = 0;
                        int compared = 0;
                        for (size_t f = 0; f < p_use; f++) {
                            AlleleType a = dataset.first[f][idx];
                            AlleleType b = dataset.first[f][jdx];
                            if (a == MISSING || b == MISSING) {
                                continue;
                            }
                            dist += std::abs(to_int(a) - to_int(b));
                            compared++;
                        }
                        if (compared == 0) {
                            continue;
                        }
                        double lbl = (double)to_int(dataset.second[jdx]);
                        if (best.size() < K) {
                            best.emplace_back(dist, lbl);
                        } else {
                            size_t worst = 0;
                            for (size_t t = 1; t < best.size(); t++) {
                                if (best[t].first > best[worst].first) {
                                    worst = t;
                                }
                            }
                            if (dist < best[worst].first) {
                                best[worst] = {dist, lbl};
                            }
                        }
                    }

                    if (best.empty()) {
                        continue;
                    }

                    double wsum = 0.0;
                    double ysum = 0.0;
                    for (const auto& nn : best) {
                        double w = 1.0 / (eps_w + (double)nn.first);
                        wsum += w;
                        ysum += w * nn.second;
                    }
                    if (wsum <= 0.0) {
                        continue;
                    }
                    double pred = ysum / wsum;
                    double err = pred - truth;
                    acc.sum_abs += std::abs(err);
                    acc.sum_sq += err * err;
                    int rounded = clamp_round_to_genotype(pred);
                    if ((double)rounded == truth) {
                        acc.correct++;
                    }
                    acc.count++;
                }
                return acc;
            };

            EvalAccum knn_accum;
            if (eval_n >= 128 && pool_threads >= 2) {
                const size_t target_tasks = std::min<std::size_t>(pool_threads * 4, 32);
                const size_t chunk = std::max<std::size_t>((eval_n + target_tasks - 1) / target_tasks, 16);
                std::vector<std::future<EvalAccum>> futures;
                for (size_t begin = 0; begin < eval_n; begin += chunk) {
                    size_t end = std::min(begin + chunk, eval_n);
                    futures.push_back(thread_pool.push([&eval_knn_range, begin, end]() { return eval_knn_range(begin, end); }));
                }
                for (auto& f : futures) {
                    auto a = f.get();
                    knn_accum.sum_abs += a.sum_abs;
                    knn_accum.sum_sq += a.sum_sq;
                    knn_accum.correct += a.correct;
                    knn_accum.count += a.count;
                }
            } else {
                knn_accum = eval_knn_range(0, eval_n);
            }

            double sum_abs_k = knn_accum.sum_abs;
            double sum_sq_k = knn_accum.sum_sq;
            size_t correct_k = knn_accum.correct;
            size_t used_k = knn_accum.count;

            if (used_k > 0) {
                knn_mae = sum_abs_k / (double)used_k;
                knn_rmse = std::sqrt(sum_sq_k / (double)used_k);
                knn_acc = (double)correct_k / (double)used_k;
            }
        }

        ImputationLooRow row;
        row.variant = (std::string)variant;
        row.n_observed = n_observed;
        row.n_missing = n_missing;
        if (rf_acc.count == 0) {
            row.oob_mae = std::numeric_limits<double>::quiet_NaN();
            row.oob_rmse = std::numeric_limits<double>::quiet_NaN();
            row.rounded_acc = std::numeric_limits<double>::quiet_NaN();
        } else {
            row.oob_mae = rf_acc.sum_abs / (double)rf_acc.count;
            row.oob_rmse = std::sqrt(rf_acc.sum_sq / (double)rf_acc.count);
            row.rounded_acc = (double)rf_acc.correct / (double)rf_acc.count;
        }

        row.ridge_loo_mae = ridge_mae;
        row.ridge_loo_rmse = ridge_rmse;
        row.ridge_rounded_acc = ridge_acc;
        row.knn_loo_mae = knn_mae;
        row.knn_loo_rmse = knn_rmse;
        row.knn_rounded_acc = knn_acc;
        loo_rows.push_back(std::move(row));
    }

    TreeBuilder PredictingHandler::make_tree_builder(const std::pair<Features, Labels>& dataset) {
        size_t mtry = ceil(sqrt(dataset.first.size()));
        return {dataset.first, dataset.second, mtry};
    }

    Window::Window(size_t max_size) :max_size(max_size), start(0) {}

    void Window::clear() {
        features.clear();
        variants.clear();
        start = 0;
    }

    std::pair<Features, Labels> Window::dataset(const Variant& v) {
        Features fs;
        Labels lbls;
        size_t none = std::numeric_limits<size_t>::max();
        size_t curr_num = none;
        bool flipped = false;

        for (size_t i = 0; i < variants.size(); i++) {
            if (variants[i] == v || variants[i] == v.reversed()) {
                curr_num = i;
                if (variants[i] == v.reversed()) {
                    flipped = true;
                }
            }
        }
        if (curr_num == none) {
            return {};
        }
        lbls = features[curr_num]->vector(flipped);
        for (size_t i = 0; i < features.size(); i++) {
            if (i != curr_num) {
                try {
                    std::vector<AlleleType> row = features[i]->vector(flipped);
                    fs.push_back(std::move(row));
                } catch (ParserException& e) {}
            }
        }
        return {std::move(fs), std::move(lbls)};
    }

    void Window::add(std::shared_ptr<AlleleVector>& alleles, const Variant& variant) {
        if (features.size() < max_size) {
            variants.push_back(variant);
            features.push_back(alleles);
        } else {
            variants.pop_front();
            features.pop_front();
            variants.push_back(variant);
            features.push_back(alleles);
        }
    }

    int Window::middle_point() {
        if (variants.empty()) {
            return -1;
        }
        return variants[variants.size() / 2].position().position();
    }

    bool Window::is_full() {
        return features.size() == max_size;
    }
}
