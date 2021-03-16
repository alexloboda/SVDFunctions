#' @useDynLib SVDFunctions
#' @importFrom Rcpp sourceCpp
NULL

#' Perform quality control on a set of allele counts
#' 
#' The methods checks which sets of allele counts pass standard quality 
#' control filters(minor allele frequency, minor allele count, HWE)
#' @param countsMatrix matrix 3-column integer matrix
#' @param maf numeric minor allele frequency
#' @param mac integer minor allele count
#' @param significance numeric significance level for HWE
#' @return logical vector(each value - whether or not the corresponding set 
#' of allele counts passed the quality control)
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

#' Select a set of controls that matches to a set of cases
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
#' @param meanGenotype Mean genotype vector
#' @param caseCounts Matrix with summary genotype counts from cases
#' @param controlsClustering cluster names for controls
#' @param minLambda Minimum possible lambda
#' @param softMinLambda Desirable minimum for lambda
#' @param softMaxLambda Desirable maximum for lambda
#' @param maxLambda Maximum possible lambda
#' @param min Minimal size of a control set that is permitted for return
#' @export
selectControls <- function (genotypeMatrix, originalGenotypeMatrix, SVDReference, 
                            meanGenotype, caseCounts, controlsClustering = NULL, minLambda = 0.75, 
                            softMinLambda = 0.9, softMaxLambda = 1.05, maxLambda = 1.3, 
                            min = 500, SV = 3) {
  stopifnot(is.matrix(genotypeMatrix))
  stopifnot(is.matrix(originalGenotypeMatrix))
  stopifnot(dim(genotypeMatrix) == dim(originalGenotypeMatrix))
  mode(genotypeMatrix) <- "numeric"
  mode(originalGenotypeMatrix) <- "integer"
  stopifnot(all(!is.na(genotypeMatrix)))
  if (nrow(genotypeMatrix) != nrow(caseCounts) || nrow(genotypeMatrix) != 
      nrow(SVDReference)) {
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
  contMeans <- rowMeans(genotypeMatrix)
  caseMeans <- apply(caseCounts, 1, function(x){
    (x[2] + 2 * x[3]) / sum(x)
  })
  
  caseCounts <- as.matrix(caseCounts)
  gmatrix <- originalGenotypeMatrix
  result <- select_controls_cpp(gmatrix, 
                                genotypeMatrix, 
                                meanGenotype, 
                                SVDReference, 
                                caseCounts, 
                                cl, 
                                stats::qchisq(stats::ppoints(1e+05), df = 1), 
                                minLambda, 
                                softMinLambda, maxLambda, softMaxLambda, min)
  result$residuals <- setNames(residuals, colnames(gmatrix))
  if (result$controls > 0) {
    result$controls <- colnames(gmatrix)[sel]
  }
  else {
    result$controls <- c()
  }
  result
}
