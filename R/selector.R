#' @useDynLib SVDFunctions
#' @importFrom Rcpp sourceCpp
NULL

#' Selection of the optimal set of controls
#' 
#' Calls C++ routine to find optimal set of controls satisfying 
#' \eqn{\lambda_GC<1.05} or if none exists – will select a set of controls with
#'  minimal \eqn{\lambda} GC satisfying \eqn{\lambda < 1.3}.
#' Otherwise no results will be returned. Minimal size of control set 
#' is \code{min} samples for privacy preservation reasons.
#' @param gmatrix Genotype matrix
#' @param svdReference Reference basis of the left singular vectors
#' @param nSV Number of singular vectors to be used for reconstruction of the 
#' @param caseCounts Matrix with summary genotype counts from cases
#' @param min Minimal size of a control set that is permitted for return
#' @param binSize sliding window size for optimal lambda search
#' @export
SelectControls <- function(gmatrix, svdReference, caseCounts, 
                            min = 500, nSV = 5, binSize = 1) {
  residuals <- ParallelResidEstimate(gmatrix, svdReference, nSV)
  control_names <- names(residuals)[order(residuals)] 
  if (dim(gmatrix)[1] != dim(caseCounts)[1] | length(residuals) != dim(gmatrix)[2]) {
    stop("Check dimensions of the matrices")
  }
  gmatrix <- as.matrix(gmatrix)
  residuals <- as.numeric(residuals)
  caseCounts <- as.matrix(caseCounts)
  result <- select_controls_cpp(gmatrix, residuals, caseCounts, 
                      stats::qchisq(ppoints(100000), df = 1), min, binSize)
  if (result$controls >= 1) {
    result$controls <- control_names[1:result$controls]
  } else {
    result$controls <- c()
  }
  result
}
