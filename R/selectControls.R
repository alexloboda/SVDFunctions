#' @useDynLib SVDFunctions
#' @importFrom Rcpp sourceCpp
NULL

#' Perform quality control on a set of allele counts.
#' 
#' The methods checks which sets of allele counts pass standard quality 
#' control filters(minor allele frequency, minor allele count, HWE).
#' @param countsMatrix matrix 3-column integer matrix.
#' @param maf numeric minor allele frequency.
#' @param mac integer minor allele count.
#' @param significance numeric significance level for HWE.
#' @return logical vector(each value - whether or not the corresponding set 
#' of allele counts passed the quality control).
#' @export
checkAlleleCounts <- function(countsMatrix, maf = 0.05, mac = 10, 
                              significance = 1e-4) {
  maf <- as.numeric(maf)
  mac <- as.integer(mac)
  significance <- as.numeric(significance)
  
  stopifnot(ncol(countsMatrix) == 3) 
  countsMatrix <- matrix(as.integer(countsMatrix), ncol = 3)
  stopifnot(all(!is.na(countsMatrix)))
  chisq_threshold <- stats::qchisq(1 - significance, 1)
  quality_control_impl(countsMatrix, maf, mac, chisq_threshold)
}

#' Select a set of controls that matches to a set of cases.
#' 
#' Finds an optimal set of controls satisfying 
#' \eqn{\lambda_GC < softmax_lambda} and \eqn{\lambda_GC > softmin_lambda}
#' or if none exists – will select a set of controls with
#'  closest \eqn{\lambda_GC} to the range \eqn{[softmin_lambda; softmax_lambda]} 
#'  satisfying \eqn{\lambda_GC < max_lambda} and \eqn{\lambda_GC > min_lambda}.
#' Otherwise no results will be returned. Minimal size of control set 
#' is \code{min} samples for privacy preservation reasons.
#' @param genotypeMatrix numeric matrix where rows are variants and columns are
#' samples. The missing values should be imputed prior calling this function.
#' @param originalGenotypeMatrix genotype matrix with missing values.
#' @param casesPDs numeric matrix where rows are variants and columns are 
#' principal directions of the data in case dataset.
#' @param casesMean numeric vector representing mean per-variant genotype value.
#' @param SVDReference reference basis of the left singular vectors.
#' @param controlsMean mean value of the reference genotypes.
#' @param caseCounts matrix with summary genotype counts from cases.
#' @param controlsClustering cluster names for controls.
#' @param minLambda minimum possible lambda.
#' @param softMinLambda desirable minimum for lambda.
#' @param softMaxLambda desirable maximum for lambda.
#' @param maxLambda maximum possible lambda.
#' @param min minimal number of clusters from control dataset to be returned.
#' @param max maximum number of clusters from control dataset to be returned.
#' @param step perform matching with the step.
#' @param iterations number of simulated annealing iterations per each subset
#' @param minCallRate numeric minimal call rate for SNP to be considered.   
#' @param method string approximation method for SA statistic, either
#'   \code{"exact"}, \code{"nystrom"}, or \code{"hybrid"}. In hybrid mode,
#'   a fixed-size Nystrom approximation is combined with a precomputed RFF ladder
#'   of increasing prefix sizes; a single SA acceptance draw is first tested
#'   against the calibrated interval from the primary Nystrom estimate and then,
#'   only if needed, against progressively larger auxiliary RFF prefixes until
#'   the sampled SA accept/reject decision is resolved within the requested
#'   CI-width threshold. A small fraction of cheap-resolved swaps is still
#'   audited exactly to keep the calibration unbiased, and an additional shadow
#'   audit stream can record full-ladder plus exact diagnostics without changing
#'   the actual SA trajectory.
#' @param n_features integer number of Nyström features (landmarks) used when
#'   \code{method = "nystrom"} or \code{method = "hybrid"}.
#' @param ciWidthThreshold numeric maximum acceptable CI width for a cheap
#'   acceptance-probability estimate. In hybrid mode, each cheap estimator uses
#'   its own exact-vs-cheap delta residuals to form a conservative acceptance
#'   interval, and exact evaluation is skipped only when some cheap estimator
#'   resolves the sampled SA accept/reject decision with CI width not exceeding
#'   this value.
#' @param calibrationQuantile numeric quantile of the exact-vs-cheap absolute
#'   delta error used to calibrate each estimator's CI half-width in hybrid mode.
#' @param minCalibrationSamples integer minimum number of successful exact
#'   audits required before using empirical calibration residuals.
#' @param maxCalibrationHistory integer maximum number of exact-audit residuals
#'   retained per estimator in hybrid mode.
#' @return A list with selected controls, lambda summaries, and \code{sa_diagnostics},
#'   which reports per-target-subset SA counters such as total swaps,
#'   CI-resolved swaps, how many of those were resolved by the primary Nyström
#'   estimate versus each auxiliary RFF ladder level, exact fallback evaluations,
#'   calibration sample counts, mean CI widths, the mean number of auxiliary
#'   ladder levels scanned before early stop, shadow-audit mismatch/error
#'   summaries, and per-temperature-bin swap statistics.
#' @export
selectControls <- function (genotypeMatrix, originalGenotypeMatrix, casesPDs, 
                            casesMean, SVDReference, controlsMean, caseCounts, 
                            controlsClustering = NULL, minLambda = 0.75, 
                            softMinLambda = 0.9, softMaxLambda = 1.05, maxLambda = 1.3, 
                            min = 500, max = 1000, step = 50, iterations = 100000, 
                            minCallRate = 0.98,
                            method = c("exact", "nystrom", "hybrid"),
                            n_features = 1024L,
                            ciWidthThreshold = 0.5,
                            calibrationQuantile = 0.9,
                            minCalibrationSamples = 8L,
                            maxCalibrationHistory = 256L) {
  iterations <- as.integer(iterations)
  method <- match.arg(method)
  n_features <- as.integer(n_features)
  ciWidthThreshold <- as.numeric(ciWidthThreshold)
  calibrationQuantile <- as.numeric(calibrationQuantile)
  minCalibrationSamples <- as.integer(minCalibrationSamples)
  maxCalibrationHistory <- as.integer(maxCalibrationHistory)
  stopifnot(iterations > 0)
  stopifnot(n_features > 0)
  stopifnot(ciWidthThreshold > 0, ciWidthThreshold <= 1)
  stopifnot(calibrationQuantile > 0, calibrationQuantile < 1)
  stopifnot(minCalibrationSamples > 0)
  stopifnot(maxCalibrationHistory > 0)
  stopifnot(is.matrix(genotypeMatrix))
  stopifnot(is.matrix(originalGenotypeMatrix))
  stopifnot(dim(genotypeMatrix) == dim(originalGenotypeMatrix))
  mode(genotypeMatrix) <- "numeric"
  mode(originalGenotypeMatrix) <- "integer"
  stopifnot(all(!is.na(genotypeMatrix)))
  if (nrow(genotypeMatrix) != nrow(caseCounts)) {
    stop("Check dimensions of the matrices")
  }
  cl <- controlsClustering
  if (is.null(cl)) {
    cl <- 0:(ncol(genotypeMatrix) - 1)
  }
  else {
    cl <- as.integer(as.factor(cl)) - 1
  }
  stopifnot(all(!is.na(cl)))
  
  names(controlsMean) <- rownames(SVDReference)
  controlsMean <- controlsMean[rownames(genotypeMatrix)]
  SVDReference <- SVDReference[rownames(genotypeMatrix), ]
  transition <- pracma::pinv(SVDReference)
  rm(SVDReference)
  
  genotypeMatrix <- genotypeMatrix - controlsMean
  genotypeMatrix <- transition %*% genotypeMatrix
  
  caseCounts <- as.matrix(caseCounts)
  gmatrix <- originalGenotypeMatrix
  
  result <- select_controls_cpp(gmatrix, 
                                genotypeMatrix, 
                                casesMean, 
                                casesPDs, 
                                caseCounts, 
                                cl, 
                                stats::qchisq(stats::ppoints(1e+07), df = 1), 
                                minLambda, 
                                softMinLambda, maxLambda, softMaxLambda, min, 
                                max, step, iterations, minCallRate,
                                method, n_features,
                                ciWidthThreshold, calibrationQuantile,
                                minCalibrationSamples, maxCalibrationHistory)
  if (length(result$controls) > 0) {
    result$controls <- colnames(gmatrix)[result$controls]
  }
  else {
    result$controls <- c()
  }
  result
}
