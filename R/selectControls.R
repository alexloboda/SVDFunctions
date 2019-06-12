#' @useDynLib SVDFunctions
#' @importFrom Rcpp sourceCpp
NULL

#' Perform quality control on a set of allele counts
#' 
#' The methods checks which sets of allele counts pass standard quality 
#' control filters(minor allele frequency > 0.05, AC > 10, HWE)
#' @param countsMatrix matrix 3-column integer matrix
#' @return logical vector(each value - whether or not the corresponding set 
#' of allele counts passed the quality control)
#' @export
checkAlleleCounts <- function(countsMatrix) {
  stopifnot(ncol(countsMatrix) == 3) 
  countsMatrix <- matrix(as.integer(countsMatrix), ncol = 3)
  stopifnot(all(!is.na(countsMatrix)))
  quality_control_impl(countsMatrix)
}

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
#' @param originalGenotypeMatrix Genotype matrix with no imputation applied 
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
selectControls <- function(genotypeMatrix, originalGenotypeMatrix,
                           SVDReference, caseCounts, 
                           minLambda = 0.75, softMinLambda = 0.9, 
                           softMaxLambda = 1.05, maxLambda = 1.3, 
                           min = 500, nSV = 5, binSize = 1) {
  stopifnot(is.matrix(genotypeMatrix))
  stopifnot(is.matrix(originalGenotypeMatrix))
  stopifnot(dim(genotypeMatrix) == dim(originalGenotypeMatrix))
  mode(genotypeMatrix) <- "numeric"
  mode(originalGenotypeMatrix) <- "numeric"
  stopifnot(any(is.na(genotypeMatrix)))
  
  if (nrow(genotypeMatrix) != nrow(caseCounts) || 
      nrow(genotypeMatrix) != nrow(SVDReference)) {
    stop("Check dimensions of the matrices")
  }
  
  residuals <- parallelResidEstimate(genotypeMatrix, SVDReference, nSV)
  new_order <- order(residuals)
  control_names <- names(residuals)[new_order] 
  
  residuals <- as.numeric(residuals)
  caseCounts <- as.matrix(caseCounts)
  gmatrix <- originalGenotypeMatrix
  result <- select_controls_cpp(gmatrix, residuals, caseCounts, 
                      stats::qchisq(stats::ppoints(100000), df = 1), 
                      minLambda, softMinLambda, maxLambda, softMaxLambda, 
                      min, binSize)
  result$residuals <- setNames(residuals, colnames(gmatrix))
  result$residuals <- result$residuals[new_order]
  if (result$controls > 0) {
    result$controls <- control_names[1:result$controls]
  } else {
    result$controls <- c()
  }
  result
}
