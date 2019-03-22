#' @useDynLib SVDFunctions
#' @importFrom Rcpp sourceCpp
NULL

#' Selection of the optimal set of controls
#' 
#' Finds an optimal set of controls satisfying 
#' \eqn{\lambda_GC < softmax_lambda} and \eqn{\lambda_GC > softmin_lambda}
#' or if none exists – will select a set of controls with
#'  closest \eqn{\lambda_GC} to the range \eqn{[softmin_lambda; softmax_lambda]} 
#'  satisfying \eqn{\lambda_GC < max_lambda} and \eqn{\lambda_GC > min_lambda}.
#' Otherwise no results will be returned. Minimal size of control set 
#' is \code{min} samples for privacy preservation reasons.
#' @param genotypeMatrix Genotype matrix
#' @param SVDReference Reference basis of the left singular vectors
#' @param caseCounts Matrix with summary genotype counts from cases
#' @param minLambda Minimum possible lambda
#' @param softMinLambda Desirable minimum for lambda
#' @param softMaxLambda Desirable maximum for lambda
#' @param maxLambda Maximum possible lambda
#' @param nSV Number of singular vectors to be used for reconstruction of the 
#' @param min Minimal size of a control set that is permitted for return
#' @param binSize sliding window size for optimal lambda search
#' @export
selectControls <- function(genotypeMatrix, SVDReference, caseCounts, 
                           minLambda = 0.75, softMinLambda = 0.9, 
                           softMaxLambda = 1.05, maxLambda = 1.3, 
                           min = 500, nSV = 5, binSize = 1) {
  gmatrix <- genotypeMatrix
  residuals <- parallelResidEstimate(gmatrix, SVDReference, nSV)
  new_order <- order(residuals)
  control_names <- names(residuals)[new_order] 
  if (dim(gmatrix)[1] != dim(caseCounts)[1] | length(residuals) != dim(gmatrix)[2]) {
    stop("Check dimensions of the matrices")
  }
  gmatrix <- as.matrix(gmatrix)
  residuals <- as.numeric(residuals)
  caseCounts <- as.matrix(caseCounts)
  result <- select_controls_cpp(gmatrix, residuals, caseCounts, 
                      stats::qchisq(stats::ppoints(100000), df = 1), 
                      minLambda, softMinLambda, maxLambda, softMaxLambda, 
                      min, binSize)
  result$residuals <- setNames(residuals, colnames(gmatrix))
  result$residuals <- result$residuals[new_order]
  if (result$controls >= 1) {
    result$controls <- control_names[1:result$controls]
  } else {
    result$controls <- c()
  }
  result
}
