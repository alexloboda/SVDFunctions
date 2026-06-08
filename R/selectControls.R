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
#' @param originalGenotypeMatrix integer genotype matrix with missing values.
#' The matrix must already use integer storage mode; \\code{selectControls}
#' throws an error instead of coercing it to avoid an additional full copy.
#' @param casesPDs numeric matrix where rows are variants and columns are 
#' principal directions of the data in case dataset.
#' @param casesMean numeric vector representing mean per-variant genotype value.
#' @param SVDReference reference basis of the left singular vectors.
#' @param controlsMean mean value of the reference genotypes.
#' @param caseCounts matrix with summary genotype counts from cases.
#' @param controlsClustering cluster names for controls. Each cluster must
#' contain at most 255 samples because per-cluster allele counts are stored
#' in one byte during matching.
#' @param minLambda minimum possible lambda.
#' @param softMinLambda desirable minimum for lambda.
#' @param softMaxLambda desirable maximum for lambda.
#' @param maxLambda maximum possible lambda.
#' @param min minimal number of clusters from control dataset to be returned.
#' @param max maximum number of clusters from control dataset to be returned.
#' @param step perform matching with the step.
#' @param iterations number of simulated annealing iterations per each subset
#' @param minCallRate numeric minimal call rate for SNP to be considered.   
#' @param saThreads optional integer number of threads for the simulated
#' annealing restarts. By default uses all available hardware threads.
#' @param exactPrecomputeThreads optional integer number of threads used to
#' precompute exact cluster-level Mahalanobis aggregates. By default uses all
#' available hardware threads.
#' @param exactClusterTileSize optional integer tile size used for exact
#' blocked aggregation across samples inside each cluster pair.
#' size.
#' @export
selectControls <- function (genotypeMatrix, originalGenotypeMatrix, casesPDs, 
                            casesMean, SVDReference, controlsMean, caseCounts, 
                            controlsClustering = NULL, minLambda = 0.75, 
                            softMinLambda = 0.9, softMaxLambda = 1.05, maxLambda = 1.3, 
                            min = 500, max = 1000, step = 50, iterations = 100000, 
                            minCallRate = 0.98, saThreads = NULL,
                            exactPrecomputeThreads = NULL,
                            exactClusterTileSize = 32L) {
  iterations <- as.integer(iterations)
  stopifnot(iterations > 0)
  saThreads <- if (is.null(saThreads)) 0L else as.integer(saThreads)
  exactPrecomputeThreads <- if (is.null(exactPrecomputeThreads)) 0L else as.integer(exactPrecomputeThreads)
  exactClusterTileSize <- as.integer(exactClusterTileSize)
  stopifnot(saThreads >= 0L)
  stopifnot(exactPrecomputeThreads >= 0L)
  stopifnot(exactClusterTileSize > 0L)
  stopifnot(is.matrix(genotypeMatrix))
  stopifnot(is.matrix(originalGenotypeMatrix))
  stopifnot(dim(genotypeMatrix) == dim(originalGenotypeMatrix))
  mode(genotypeMatrix) <- "numeric"
  if (!is.integer(originalGenotypeMatrix)) {
    stop("originalGenotypeMatrix must already be stored as an integer matrix")
  }
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
  if (any(tabulate(cl + 1L) > 255L)) {
    stop("Each control cluster must contain at most 255 samples")
  }
  
  names(controlsMean) <- rownames(SVDReference)
  controlsMean <- controlsMean[rownames(genotypeMatrix)]
  SVDReference <- SVDReference[rownames(genotypeMatrix), ]
  transition <- pracma::pinv(SVDReference)
  rm(SVDReference)
  
  meanOffset <- as.vector(transition %*% controlsMean)
  genotypeMatrix <- transition %*% genotypeMatrix
  genotypeMatrix <- genotypeMatrix - meanOffset
  
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
                                saThreads, exactPrecomputeThreads,
                                exactClusterTileSize)
  if (length(result$controls) > 0) {
    result$controls <- colnames(gmatrix)[result$controls]
  }
  else {
    result$controls <- c()
  }
  result
}
