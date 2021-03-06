# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

quality_control_impl <- function(case_counts, maf, mac, chi2boundary) {
    .Call('_SVDFunctions_quality_control_impl', PACKAGE = 'SVDFunctions', case_counts, maf, mac, chi2boundary)
}

select_controls_cpp <- function(gmatrix, residuals, cc, clustering, chi2fn, min_lambda, lb_lambda, max_lambda, ub_lambda, min) {
    .Call('_SVDFunctions_select_controls_cpp', PACKAGE = 'SVDFunctions', gmatrix, residuals, cc, clustering, chi2fn, min_lambda, lb_lambda, max_lambda, ub_lambda, min)
}

parse_vcf <- function(filename, samples, bad_positions, variants, DP, GQ, gmatrix, predictMissing, regions, binary_prefix, missingRateThreshold) {
    .Call('_SVDFunctions_parse_vcf', PACKAGE = 'SVDFunctions', filename, samples, bad_positions, variants, DP, GQ, gmatrix, predictMissing, regions, binary_prefix, missingRateThreshold)
}

parse_binary_file <- function(variants, samples, regions, binary_file, metafile, r_min_maf, r_max_maf, r_min_cr, r_min_mac, r_max_mac, report_singletons, requiredDP, requiredGQ) {
    .Call('_SVDFunctions_parse_binary_file', PACKAGE = 'SVDFunctions', variants, samples, regions, binary_file, metafile, r_min_maf, r_max_maf, r_min_cr, r_min_mac, r_max_mac, report_singletons, requiredDP, requiredGQ)
}

