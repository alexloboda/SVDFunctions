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
#' @param residuals Vector of relative residual vector norms for all prospective controls
#' @param case_counts Matrix with summary genotype counts from cases
#' @param min Minimal size of a control set that is permitted for return
#' @export
select_controls <- function(gmatrix, residuals, case_counts, 
                            min = 500) {
  if (dim(gmatrix)[1] != dim(case_counts)[1] | length(residuals) != dim(gmatrix)[2]) {
    stop("Check dimensions of the matrices")
  }
  gmatrix <- as.matrix(gmatrix)
  residuals <- as.numeric(residuals)
  case_counts <- as.matrix(case_counts)
  out <- select_controls_cpp(gmatrix, residuals, case_counts, 
                      stats::qchisq(ppoints(100000), df = 1), min)
  out$lambda <- setNames(out$lambda, as.character(min:(min + length(out$lambda) - 1)))
  out
}