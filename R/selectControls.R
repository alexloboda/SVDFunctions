#' @useDynLib SVDFunctions, .registration = TRUE
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

prepareControlsClustering <- function(controlsClustering, controlNames,
                                      nControls = length(controlNames)) {
  if (is.null(controlsClustering)) {
    if (is.null(controlNames)) {
      stop(paste("genotypeMatrix must have column names when",
                 "controlsClustering is NULL"))
    }
    return(list(sampleClusterIds = seq_len(nControls) - 1L,
                clusterLabels = controlNames))
  }
  if (length(controlsClustering) != nControls) {
    stop("controlsClustering must have one entry per control sample")
  }

  factoredClustering <- as.factor(controlsClustering)
  list(sampleClusterIds = as.integer(factoredClustering) - 1L,
       clusterLabels = levels(factoredClustering))
}

resolveSelectedControls <- function(clusterIds, sampleClusterIds, clusterLabels,
                                    controlNames, returnClusters = FALSE) {
  if (isTRUE(returnClusters)) {
    return(clusterLabels[clusterIds])
  }
  if (is.null(controlNames)) {
    return(clusterLabels[clusterIds])
  }

  samplesByCluster <- split(controlNames, sampleClusterIds)
  unlist(samplesByCluster[as.character(clusterIds - 1L)], use.names = FALSE)
}

#' Select a set of controls that matches to a set of cases.
#' 
#' Finds an optimal set of controls satisfying 
#' \eqn{\lambda_GC < softmax_lambda} and \eqn{\lambda_GC > softmin_lambda}
#' or if none exists -- will select a set of controls with
#'  closest \eqn{\lambda_GC} to the range \eqn{[softmin_lambda; softmax_lambda]} 
#'  satisfying \eqn{\lambda_GC < max_lambda} and \eqn{\lambda_GC > min_lambda}.
#' Otherwise no results will be returned. Minimal size of control set 
#' is \code{min} samples for privacy preservation reasons.
#' @param genotypeMatrix numeric matrix where rows are variants and columns are
#' samples. The missing values should be imputed prior calling this function.
#' Rows are identified by their variant names (\code{rownames}) and may contain
#' a superset of \code{caseVariants}.
#' @param originalGenotypeMatrix integer genotype matrix with missing values.
#' The matrix must already use integer storage mode; \code{selectControls}
#' throws an error instead of coercing it to avoid an additional full copy. Its
#' rows must correspond to those of \code{genotypeMatrix}.
#' @param casesPDs numeric matrix of case principal directions in the PCA-like
#' space (rows are reduced components, columns are principal directions).
#' @param casesMean numeric vector of case means in the PCA-like space, one
#' value per reduced component.
#' @param SVDReference reference basis of the left singular vectors.
#' @param controlsMean mean value of the reference genotypes.
#' @param caseCounts matrix with summary genotype counts from cases, one row per
#' case variant in the order given by \code{caseVariants}.
#' @param caseVariants optional character vector naming the case variants in the
#' row order of \code{caseCounts}. Every entry must occur in the row names of
#' \code{genotypeMatrix} and \code{SVDReference}. Defaults to the control variant
#' names, i.e. cases and controls share the same variants.
#' @param controlsClustering cluster names for controls, one per matrix column.
#' Each cluster must contain at most 255 samples because per-cluster allele
#' counts are stored in one byte during matching. This argument is required
#' when \code{genotypeMatrix} has no column names.
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
#' @param returnClusters logical; controls the form of the returned
#' \code{controls} element. When \code{FALSE} (default) it contains sample
#' names. When \code{TRUE} it contains the
#' cluster identifiers instead -- the original labels supplied via
#' \code{controlsClustering}, or sample names when no clusters were provided.
#' Cluster identifiers are also returned when the controls have no column names.
#' @return a list with the matching diagnostics (\code{lambda},
#' \code{optimal_lambda}, \code{statistics}, \code{pvals}, \code{snps}) and a
#' \code{controls} element holding either the selected sample names or, when
#' \code{returnClusters = TRUE}, the selected cluster identifiers.
#' @export
selectControls <- function (genotypeMatrix, originalGenotypeMatrix, casesPDs, 
                            casesMean, SVDReference, controlsMean, caseCounts, 
                            caseVariants = NULL,
                            controlsClustering = NULL, minLambda = 0.75, 
                            softMinLambda = 0.9, softMaxLambda = 1.05, maxLambda = 1.3, 
                            min = 500, max = 1000, step = 50, iterations = 100000, 
                            minCallRate = 0.98, saThreads = NULL,
                            exactPrecomputeThreads = NULL,
                            exactClusterTileSize = 32L,
                            returnClusters = FALSE) {
  iterations <- as.integer(iterations)
  stopifnot(iterations > 0)
  returnClusters <- isTRUE(returnClusters)
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
  controlVariants <- rownames(genotypeMatrix)
  if (is.null(controlVariants)) {
    stop("genotypeMatrix must have variant row names")
  }
  if (is.null(caseVariants)) {
    caseVariants <- controlVariants
  }
  caseVariants <- as.character(caseVariants)
  if (anyDuplicated(caseVariants)) {
    stop("caseVariants must not contain duplicates")
  }
  if (nrow(caseCounts) != length(caseVariants)) {
    stop("caseCounts must have one row per case variant")
  }
  variantRows <- match(caseVariants, controlVariants)
  if (anyNA(variantRows)) {
    stop("caseVariants must be a subset of the control genotype matrix rows")
  }
  clusteringInfo <- prepareControlsClustering(controlsClustering,
                                              colnames(genotypeMatrix),
                                              ncol(genotypeMatrix))
  cl <- clusteringInfo$sampleClusterIds
  clusterLabels <- clusteringInfo$clusterLabels
  stopifnot(all(!is.na(cl)))
  if (any(tabulate(cl + 1L) > 255L)) {
    stop("Each control cluster must contain at most 255 samples")
  }
  
  names(controlsMean) <- rownames(SVDReference)
  if (!all(caseVariants %in% rownames(SVDReference))) {
    stop("Every case variant must occur in SVDReference")
  }
  transition <- pinv(SVDReference[caseVariants, , drop = FALSE])
  meanOffset <- as.vector(transition %*% controlsMean[caseVariants])
  rm(SVDReference)
  
  # Project the controls into the case PCA-like space through an index view over
  # the control variants instead of slicing genotypeMatrix into a large copy:
  # the transition is widened to every control variant with zero columns for
  # those absent from the cases, so multiplying the full matrix yields the same
  # reduced result as projecting only the shared variants.
  if (length(variantRows) != length(controlVariants) ||
      any(variantRows != seq_along(controlVariants))) {
    widened <- matrix(0, nrow(transition), length(controlVariants))
    widened[, variantRows] <- transition
    transition <- widened
  }
  genotypeMatrix <- transition %*% genotypeMatrix
  genotypeMatrix <- genotypeMatrix - meanOffset
  
  caseCounts <- as.matrix(caseCounts)
  gmatrix <- originalGenotypeMatrix
  
  result <- select_controls_cpp(gmatrix, 
                                genotypeMatrix, 
                                casesMean, 
                                casesPDs, 
                                caseCounts, 
                                variantRows - 1L, 
                                cl, 
                                stats::qchisq(stats::ppoints(1e+07), df = 1), 
                                minLambda, 
                                softMinLambda, maxLambda, softMaxLambda, min, 
                                max, step, iterations, minCallRate,
                                saThreads, exactPrecomputeThreads,
                                exactClusterTileSize)
  if (length(result$controls) > 0) {
    result$controls <- resolveSelectedControls(result$controls, cl,
                                               clusterLabels, colnames(gmatrix),
                                               returnClusters)
  }
  else {
    result$controls <- c()
  }
  result
}
