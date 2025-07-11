// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// quality_control_impl
LogicalVector quality_control_impl(const IntegerMatrix& case_counts, const NumericVector& maf, const IntegerVector& mac, const NumericVector& chi2boundary);
RcppExport SEXP _SVDFunctions_quality_control_impl(SEXP case_countsSEXP, SEXP mafSEXP, SEXP macSEXP, SEXP chi2boundarySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type case_counts(case_countsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type maf(mafSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type mac(macSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type chi2boundary(chi2boundarySEXP);
    rcpp_result_gen = Rcpp::wrap(quality_control_impl(case_counts, maf, mac, chi2boundary));
    return rcpp_result_gen;
END_RCPP
}
// subsample_mvn
List subsample_mvn(NumericMatrix& matrix, IntegerVector size, NumericVector& mean, NumericMatrix& cov);
RcppExport SEXP _SVDFunctions_subsample_mvn(SEXP matrixSEXP, SEXP sizeSEXP, SEXP meanSEXP, SEXP covSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type matrix(matrixSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type cov(covSEXP);
    rcpp_result_gen = Rcpp::wrap(subsample_mvn(matrix, size, mean, cov));
    return rcpp_result_gen;
END_RCPP
}
// select_controls_cpp
List select_controls_cpp(IntegerMatrix& gmatrix, NumericMatrix& gmatrix_rs, NumericVector& mean, NumericMatrix& directions, IntegerMatrix& cc, IntegerVector& clustering, NumericVector& chi2fn, double min_lambda, double lb_lambda, double max_lambda, double ub_lambda, int min, int max, int step, int sa_iterations, double min_call_rate);
RcppExport SEXP _SVDFunctions_select_controls_cpp(SEXP gmatrixSEXP, SEXP gmatrix_rsSEXP, SEXP meanSEXP, SEXP directionsSEXP, SEXP ccSEXP, SEXP clusteringSEXP, SEXP chi2fnSEXP, SEXP min_lambdaSEXP, SEXP lb_lambdaSEXP, SEXP max_lambdaSEXP, SEXP ub_lambdaSEXP, SEXP minSEXP, SEXP maxSEXP, SEXP stepSEXP, SEXP sa_iterationsSEXP, SEXP min_call_rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix& >::type gmatrix(gmatrixSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type gmatrix_rs(gmatrix_rsSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type directions(directionsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix& >::type cc(ccSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type clustering(clusteringSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type chi2fn(chi2fnSEXP);
    Rcpp::traits::input_parameter< double >::type min_lambda(min_lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type lb_lambda(lb_lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type max_lambda(max_lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type ub_lambda(ub_lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type min(minSEXP);
    Rcpp::traits::input_parameter< int >::type max(maxSEXP);
    Rcpp::traits::input_parameter< int >::type step(stepSEXP);
    Rcpp::traits::input_parameter< int >::type sa_iterations(sa_iterationsSEXP);
    Rcpp::traits::input_parameter< double >::type min_call_rate(min_call_rateSEXP);
    rcpp_result_gen = Rcpp::wrap(select_controls_cpp(gmatrix, gmatrix_rs, mean, directions, cc, clustering, chi2fn, min_lambda, lb_lambda, max_lambda, ub_lambda, min, max, step, sa_iterations, min_call_rate));
    return rcpp_result_gen;
END_RCPP
}
// parse_vcf
List parse_vcf(const CharacterVector& filename, const CharacterVector& samples, const CharacterVector& bad_positions, const CharacterVector& variants, const IntegerVector& DP, const IntegerVector& GQ, const LogicalVector& gmatrix, const LogicalVector& predictMissing, const CharacterVector& regions, const CharacterVector& binary_prefix, const NumericVector& missingRateThreshold, Rcpp::Nullable<int> seed);
RcppExport SEXP _SVDFunctions_parse_vcf(SEXP filenameSEXP, SEXP samplesSEXP, SEXP bad_positionsSEXP, SEXP variantsSEXP, SEXP DPSEXP, SEXP GQSEXP, SEXP gmatrixSEXP, SEXP predictMissingSEXP, SEXP regionsSEXP, SEXP binary_prefixSEXP, SEXP missingRateThresholdSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const CharacterVector& >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type bad_positions(bad_positionsSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type variants(variantsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type DP(DPSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type GQ(GQSEXP);
    Rcpp::traits::input_parameter< const LogicalVector& >::type gmatrix(gmatrixSEXP);
    Rcpp::traits::input_parameter< const LogicalVector& >::type predictMissing(predictMissingSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type regions(regionsSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type binary_prefix(binary_prefixSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type missingRateThreshold(missingRateThresholdSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(parse_vcf(filename, samples, bad_positions, variants, DP, GQ, gmatrix, predictMissing, regions, binary_prefix, missingRateThreshold, seed));
    return rcpp_result_gen;
END_RCPP
}
// parse_binary_file
List parse_binary_file(const CharacterVector& variants, const CharacterVector& samples, const CharacterVector& regions, const CharacterVector& binary_file, const CharacterVector& metafile, const NumericVector& r_min_maf, const NumericVector& r_max_maf, const NumericVector& r_min_cr, const IntegerVector& r_min_mac, const IntegerVector& r_max_mac, const LogicalVector& report_singletons, const IntegerVector& requiredDP, const IntegerVector requiredGQ);
RcppExport SEXP _SVDFunctions_parse_binary_file(SEXP variantsSEXP, SEXP samplesSEXP, SEXP regionsSEXP, SEXP binary_fileSEXP, SEXP metafileSEXP, SEXP r_min_mafSEXP, SEXP r_max_mafSEXP, SEXP r_min_crSEXP, SEXP r_min_macSEXP, SEXP r_max_macSEXP, SEXP report_singletonsSEXP, SEXP requiredDPSEXP, SEXP requiredGQSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const CharacterVector& >::type variants(variantsSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type regions(regionsSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type binary_file(binary_fileSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type metafile(metafileSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type r_min_maf(r_min_mafSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type r_max_maf(r_max_mafSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type r_min_cr(r_min_crSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type r_min_mac(r_min_macSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type r_max_mac(r_max_macSEXP);
    Rcpp::traits::input_parameter< const LogicalVector& >::type report_singletons(report_singletonsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type requiredDP(requiredDPSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type requiredGQ(requiredGQSEXP);
    rcpp_result_gen = Rcpp::wrap(parse_binary_file(variants, samples, regions, binary_file, metafile, r_min_maf, r_max_maf, r_min_cr, r_min_mac, r_max_mac, report_singletons, requiredDP, requiredGQ));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SVDFunctions_quality_control_impl", (DL_FUNC) &_SVDFunctions_quality_control_impl, 4},
    {"_SVDFunctions_subsample_mvn", (DL_FUNC) &_SVDFunctions_subsample_mvn, 4},
    {"_SVDFunctions_select_controls_cpp", (DL_FUNC) &_SVDFunctions_select_controls_cpp, 16},
    {"_SVDFunctions_parse_vcf", (DL_FUNC) &_SVDFunctions_parse_vcf, 12},
    {"_SVDFunctions_parse_binary_file", (DL_FUNC) &_SVDFunctions_parse_binary_file, 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_SVDFunctions(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
