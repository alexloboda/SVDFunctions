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

# Shared core of the prepareControlsSpace methods. Everything except the
# projection itself is identical between the in-memory and the file backed path,
# so the difference is confined to `project`, which receives the k x
# n_case_variants transition, the rows the case variants occupy in the controls
# and the control means of those variants, and returns the reduced coordinates.
buildControlsSpace <- function(controlVariants, sampleNames, nSamples,
                               SVDReference, controlsMean, casesMean, casesPDs,
                               caseVariants, controlsClustering, project) {
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
  variantRows <- match(caseVariants, controlVariants)
  if (anyNA(variantRows)) {
    stop("caseVariants must be a subset of the control genotype matrix rows")
  }
  clusteringInfo <- prepareControlsClustering(controlsClustering, sampleNames,
                                              nSamples)
  cl <- clusteringInfo$sampleClusterIds
  stopifnot(all(!is.na(cl)))
  if (any(tabulate(cl + 1L) > 255L)) {
    stop("Each control cluster must contain at most 255 samples")
  }

  names(controlsMean) <- rownames(SVDReference)
  if (!all(caseVariants %in% rownames(SVDReference))) {
    stop("Every case variant must occur in SVDReference")
  }
  transition <- pinv(SVDReference[caseVariants, , drop = FALSE])
  variantMean <- as.vector(controlsMean[caseVariants])
  rm(SVDReference)

  points <- project(transition, variantRows, variantMean)
  # The in-memory projection inherits these from the multiplication; set them
  # explicitly so both paths label the columns of `points` the same way.
  colnames(points) <- sampleNames

  list(points = points,
       mean = casesMean,
       cov = tcrossprod(as.matrix(casesPDs)),
       clusterIds = cl + 1L,
       clusterLabels = clusteringInfo$clusterLabels,
       sampleNames = sampleNames,
       caseVariants = caseVariants,
       variantRows = variantRows)
}

#' Project control genotypes into the case PCA-like space.
#'
#' First stage of \code{\link{selectControls}}. It reduces the control genotypes
#' to the PCA-like coordinates the cases were reduced to, and collects the target
#' mean and covariance that describe the case distribution in that space. The
#' returned list carries everything \code{\link{mvnSubsampleClusters}} needs, so
#' the annealing can be run directly on it.
#'
#' The controls can be given either as a matrix or as the path of a file written
#' by \code{\link{writeControlsMatrix}}. The file backed method streams the
#' matrix one variant at a time and never materialises it: the projection
#' decomposes into one contribution per variant, so its peak memory is the
#' reduced result plus a handful of variant rows, whatever the size of the
#' matrix on disk.
#'
#' The two methods sum the same products in a different order, so the reduced
#' coordinates agree to rounding rather than bit for bit. The search that
#' consumes them is discrete, so the two can settle on different -- equally
#' valid -- subsets of controls. Within either method the result is exact and
#' does not depend on \code{threads}.
#' @param genotypeMatrix either the numeric control genotype matrix described in
#' \code{\link{selectControls}}, or the path of a control matrix file holding it.
#' @param threads optional integer number of threads for the file backed
#' projection. By default uses all available hardware threads. Ignored for the
#' in-memory method.
#' @param ... passed to the method.
#' @inheritParams selectControls
#' @return a list with the reduced control coordinates (\code{points}, one column
#' per control sample), the target \code{mean} and \code{cov} of the case
#' distribution, the one-based per-sample \code{clusterIds} together with their
#' \code{clusterLabels}, the control \code{sampleNames}, the resolved
#' \code{caseVariants} and the \code{variantRows} they occupy in
#' \code{genotypeMatrix}.
#' @seealso \code{\link{writeControlsMatrix}}
#' @export
prepareControlsSpace <- function(genotypeMatrix, ...) {
  UseMethod("prepareControlsSpace")
}

#' @rdname prepareControlsSpace
#' @export
prepareControlsSpace.matrix <- function(genotypeMatrix, SVDReference, controlsMean,
                                        casesMean, casesPDs, caseVariants = NULL,
                                        controlsClustering = NULL, ...) {
  if (mode(genotypeMatrix) != "numeric") {
    stop("genotypeMatrix must already be stored as a numeric matrix")
  }
  stopifnot(all(!is.na(genotypeMatrix)))
  controlVariants <- rownames(genotypeMatrix)

  project <- function(transition, variantRows, variantMean) {
    meanOffset <- as.vector(transition %*% variantMean)
    # Project the controls into the case PCA-like space through an index view
    # over the control variants instead of slicing genotypeMatrix into a large
    # copy: the transition is widened to every control variant with zero columns
    # for those absent from the cases, so multiplying the full matrix yields the
    # same reduced result as projecting only the shared variants.
    if (length(variantRows) != length(controlVariants) ||
        any(variantRows != seq_along(controlVariants))) {
      widened <- matrix(0, nrow(transition), length(controlVariants))
      widened[, variantRows] <- transition
      transition <- widened
    }
    points <- transition %*% genotypeMatrix
    points - meanOffset
  }

  buildControlsSpace(controlVariants, colnames(genotypeMatrix),
                     ncol(genotypeMatrix), SVDReference, controlsMean, casesMean,
                     casesPDs, caseVariants, controlsClustering, project)
}

#' @rdname prepareControlsSpace
#' @export
prepareControlsSpace.character <- function(genotypeMatrix, SVDReference, controlsMean,
                                           casesMean, casesPDs, caseVariants = NULL,
                                           controlsClustering = NULL, threads = NULL,
                                           ...) {
  path <- path.expand(genotypeMatrix)
  info <- readControlsMatrixInfo(path)
  if (info$dtype == "clusterCounts") {
    stop("the projection needs a genotype matrix, not collapsed cluster counts")
  }
  threads <- if (is.null(threads)) 0L else as.integer(threads)
  stopifnot(threads >= 0L)

  # The C++ side centers each variant as it streams it, which is exactly the
  # meanOffset subtraction the in-memory path performs after the fact.
  project <- function(transition, variantRows, variantMean) {
    project_controls_file_cpp(path, variantRows - 1L, transition, variantMean,
                              threads)
  }

  buildControlsSpace(info$rownames, info$colnames, info$ncol, SVDReference,
                     controlsMean, casesMean, casesPDs, caseVariants,
                     controlsClustering, project)
}

#' Select multivariate normal subsamples of a set of points.
#'
#' Second stage of \code{\link{selectControls}}, and a self-contained
#' multivariate normality search: it only sees points in some reduced space, the
#' clusters they are swapped in as units, and the normal distribution to match.
#' Nothing about genotypes enters here. Simulated annealing minimizes a BHEP-type
#' normality statistic for every target subset size in turn, so the result is one
#' candidate subset per size rather than a single answer.
#' @param points numeric matrix where columns are points and rows are their
#' components, e.g. the \code{points} element of \code{\link{prepareControlsSpace}}.
#' @param mean numeric vector, one value per component, giving the mean of the
#' distribution the subsample should follow.
#' @param cov numeric covariance matrix of that distribution, with one row and
#' column per component.
#' @param clusterIds optional one-based cluster id per column of \code{points}.
#' Points sharing a cluster are added and removed as a unit. Defaults to one
#' cluster per point.
#' @param min minimal number of clusters to select.
#' @param max maximum number of clusters to select.
#' @param step increment between consecutive target sizes.
#' @param iterations number of simulated annealing iterations per target size.
#' @param seed optional integer seed for the simulated annealing. When
#' \code{NULL} (default) it is drawn from R's generator, so \code{set.seed}
#' fixes the search; see \code{\link{selectControls}} for the exact guarantee.
#' @param threads optional integer number of threads for the annealing restarts.
#' By default uses all available hardware threads.
#' @param precomputeThreads optional integer number of threads used to precompute
#' the exact cluster-level Mahalanobis aggregates. By default uses all available
#' hardware threads.
#' @param clusterTileSize optional integer tile size used for exact blocked
#' aggregation across points inside each cluster pair.
#' @return a list with \code{clusters}, a list holding the one-based cluster ids
#' of the best subset found for each target size, and \code{statistics}, the
#' normality statistic each of those subsets reached.
#' @export
mvnSubsampleClusters <- function(points, mean, cov, clusterIds = NULL,
                                 min = 500, max = 1000, step = 50,
                                 iterations = 100000, seed = NULL,
                                 threads = NULL, precomputeThreads = NULL,
                                 clusterTileSize = 32L) {
  stopifnot(is.matrix(points))
  if (mode(points) != "numeric") {
    stop("points must already be stored as a numeric matrix")
  }
  cov <- as.matrix(cov)
  mean <- as.numeric(mean)
  if (nrow(cov) != nrow(points) || ncol(cov) != nrow(points)) {
    stop("cov must be square with one row per component of points")
  }
  if (length(mean) != nrow(points)) {
    stop("mean must have one value per component of points")
  }

  if (is.null(clusterIds)) {
    clusterIds <- seq_len(ncol(points))
  }
  clusterIds <- as.integer(clusterIds)
  if (length(clusterIds) != ncol(points)) {
    stop("clusterIds must have one entry per column of points")
  }
  if (anyNA(clusterIds) || any(clusterIds < 1L)) {
    stop("clusterIds must be one-based cluster indices")
  }

  seed <- resolveSeed(seed)
  iterations <- as.integer(iterations)
  min <- as.integer(min)
  max <- as.integer(max)
  step <- as.integer(step)
  threads <- if (is.null(threads)) 0L else as.integer(threads)
  precomputeThreads <- if (is.null(precomputeThreads)) 0L else as.integer(precomputeThreads)
  clusterTileSize <- as.integer(clusterTileSize)
  stopifnot(iterations > 0)
  # A non-positive step would never advance the target size in the C++ loop.
  stopifnot(min >= 0L, max >= 0L, step > 0L)
  stopifnot(threads >= 0L)
  stopifnot(precomputeThreads >= 0L)
  stopifnot(clusterTileSize > 0L)

  mvn_subsample_clusters_cpp(points, clusterIds - 1L, mean, cov,
                             min, max, step, iterations, seed,
                             threads, precomputeThreads, clusterTileSize)
}

# The chi-squared quantile table the lambda_GC estimate is read off.
matchingChisqTable <- function() {
  stats::qchisq(stats::ppoints(1e+07), df = 1)
}

validateCandidates <- function(candidates) {
  if (!is.list(candidates) || is.null(candidates$clusters) ||
      is.null(candidates$statistics)) {
    stop("candidates must be a list with clusters and statistics elements")
  }
  candidates
}

validateClusterIds <- function(clusterIds, nSamples) {
  clusterIds <- as.integer(clusterIds)
  if (length(clusterIds) != nSamples) {
    stop("clusterIds must have one entry per control sample")
  }
  if (anyNA(clusterIds) || any(clusterIds < 1L)) {
    stop("clusterIds must be one-based cluster indices")
  }
  clusterIds
}

#' Pick the genetically best matching subset among control candidates.
#'
#' Third stage of \code{\link{selectControls}}, and the only one that looks at
#' genotypes. For every candidate subset produced by
#' \code{\link{mvnSubsampleClusters}} it recomputes the per-variant control
#' counts, applies call rate and allele count quality control, derives per-variant
#' association p-values against the case counts and from those the genomic
#' inflation factor \eqn{\lambda_GC}. It then keeps the candidates within
#' \eqn{[minLambda; maxLambda]} and prefers the ones inside
#' \eqn{[softMinLambda; softMaxLambda]}, falling back to the one closest to that
#' range.
#'
#' The controls can be given as a matrix, as the path of a raw genotype matrix
#' written by \code{\link{writeControlsMatrix}} with \code{dtype = "u8"}, or as
#' the path of a per-cluster counts file written by
#' \code{\link{collapseControlsMatrix}}. The collapsed file is the cheapest of
#' the three: it already holds exactly the counts this stage needs, so nothing is
#' recomputed and no individual genotype is read. It also carries its own
#' clustering, so \code{clusterIds} defaults to the one it was built with.
#' @param originalGenotypeMatrix either the integer control genotype matrix
#' described in \code{\link{selectControls}}, or the path of a control matrix
#' file holding raw genotypes or collapsed per-cluster counts.
#' @param candidates list as returned by \code{\link{mvnSubsampleClusters}}, with
#' a \code{clusters} element holding one-based cluster ids per candidate and a
#' matching \code{statistics} element.
#' @param clusterIds one-based cluster id per control sample, as returned by
#' \code{\link{prepareControlsSpace}}. For a collapsed counts file it defaults to
#' the clustering stored in the file.
#' @param variantRows optional one-based rows of \code{originalGenotypeMatrix}
#' holding the case variants, in the row order of \code{caseCounts}. Defaults to
#' every row, i.e. cases and controls share the same variants in the same order.
#' @param caseVariants optional names of the case variants, in the row order of
#' \code{caseCounts}. Only for the file backed method, where it is resolved
#' against the variant names stored in the file instead of \code{variantRows}.
#' @param ... passed to the method.
#' @inheritParams selectControls
#' @return a list with the matching diagnostics (\code{lambda},
#' \code{optimal_lambda}, \code{statistics}, \code{pvals}, \code{snps}) and a
#' \code{controls} element holding the one-based cluster ids of the selected
#' candidate, empty when no candidate satisfied the hard lambda bounds.
#' @seealso \code{\link{collapseControlsMatrix}}
#' @export
matchControlCandidates <- function(originalGenotypeMatrix, ...) {
  UseMethod("matchControlCandidates")
}

#' @rdname matchControlCandidates
#' @export
matchControlCandidates.matrix <- function(originalGenotypeMatrix, caseCounts,
                                          candidates, clusterIds,
                                          variantRows = NULL, minLambda = 0.75,
                                          softMinLambda = 0.9, softMaxLambda = 1.05,
                                          maxLambda = 1.3, min = 500,
                                          minCallRate = 0.98, ...) {
  if (!is.integer(originalGenotypeMatrix)) {
    stop("originalGenotypeMatrix must already be stored as an integer matrix")
  }
  validateCandidates(candidates)
  clusterIds <- validateClusterIds(clusterIds, ncol(originalGenotypeMatrix))
  caseCounts <- as.matrix(caseCounts)
  if (is.null(variantRows)) {
    variantRows <- seq_len(nrow(originalGenotypeMatrix))
  }
  variantRows <- as.integer(variantRows)
  if (nrow(caseCounts) != length(variantRows)) {
    stop("caseCounts must have one row per case variant")
  }

  match_controls_cpp(originalGenotypeMatrix, caseCounts, variantRows - 1L,
                     clusterIds - 1L,
                     lapply(candidates$clusters, as.integer),
                     as.numeric(candidates$statistics),
                     matchingChisqTable(),
                     minLambda, softMinLambda, maxLambda, softMaxLambda,
                     min, minCallRate)
}

#' @rdname matchControlCandidates
#' @export
matchControlCandidates.character <- function(originalGenotypeMatrix, caseCounts,
                                             candidates, clusterIds = NULL,
                                             variantRows = NULL, caseVariants = NULL,
                                             minLambda = 0.75, softMinLambda = 0.9,
                                             softMaxLambda = 1.05, maxLambda = 1.3,
                                             min = 500, minCallRate = 0.98, ...) {
  path <- path.expand(originalGenotypeMatrix)
  info <- readControlsMatrixInfo(path)
  collapsed <- info$dtype == "clusterCounts"
  if (!collapsed && info$dtype != "u8") {
    stop(paste("matching needs a matrix written with dtype \"u8\" or collapsed",
               "cluster counts"))
  }
  validateCandidates(candidates)

  if (is.null(clusterIds)) {
    if (collapsed) {
      if (any(info$colweights < 1L)) {
        stop("the cluster counts file holds an empty cluster")
      }
      # A collapsed file has no samples left to count, so the clustering is
      # rebuilt from the sample sizes it stores alongside the counts.
      clusterIds <- rep(seq_len(info$ncol), info$colweights)
    } else {
      clusterIds <- seq_len(info$ncol)
    }
  }
  clusterIds <- validateClusterIds(clusterIds, sum(info$colweights))

  caseCounts <- as.matrix(caseCounts)
  if (is.null(variantRows)) {
    variantRows <- controlsMatrixVariantRows(info, caseVariants)
  }
  variantRows <- as.integer(variantRows)
  if (nrow(caseCounts) != length(variantRows)) {
    stop("caseCounts must have one row per case variant")
  }

  match_controls_file_cpp(path, caseCounts, variantRows - 1L, clusterIds - 1L,
                          lapply(candidates$clusters, as.integer),
                          as.numeric(candidates$statistics),
                          matchingChisqTable(),
                          minLambda, softMinLambda, maxLambda, softMaxLambda,
                          min, minCallRate)
}

# Stages two and three plus the resolution of the selected clusters back to
# names. `match` runs the third stage over the candidates the search produced;
# it is the only part that differs between the in-memory and file backed entry
# points.
runControlSelection <- function(space, match, sampleNames, returnClusters, min, max,
                                step, iterations, seed, saThreads,
                                exactPrecomputeThreads, exactClusterTileSize) {
  candidates <- mvnSubsampleClusters(space$points, space$mean, space$cov,
                                     space$clusterIds, min = min, max = max,
                                     step = step, iterations = iterations,
                                     seed = seed, threads = saThreads,
                                     precomputeThreads = exactPrecomputeThreads,
                                     clusterTileSize = exactClusterTileSize)
  result <- match(candidates)
  if (length(result$controls) > 0) {
    # resolveSelectedControls expects zero-based per-sample cluster ids.
    result$controls <- resolveSelectedControls(result$controls,
                                               space$clusterIds - 1L,
                                               space$clusterLabels,
                                               sampleNames, returnClusters)
  }
  else {
    result$controls <- c()
  }
  result
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
#' @param seed optional integer seed for the simulated annealing. When
#' \code{NULL} (default) it is drawn from R's generator, so \code{set.seed}
#' fixes the selection. The same seed reproduces the same controls exactly for a
#' given build of the package on a given machine, and does so independently of
#' \code{saThreads} and \code{exactPrecomputeThreads}. It is not a
#' cross-platform guarantee: rebuilding with different compiler flags or a
#' different standard library changes the result, and because the annealing is a
#' discrete search, any such difference yields an entirely different subset
#' rather than a slightly different one.
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
                            returnClusters = FALSE, seed = NULL) {
  seed <- resolveSeed(seed)
  returnClusters <- isTRUE(returnClusters)
  stopifnot(is.matrix(genotypeMatrix))
  stopifnot(is.matrix(originalGenotypeMatrix))
  stopifnot(dim(genotypeMatrix) == dim(originalGenotypeMatrix))
  # Checked here rather than left to matchControlCandidates so that a wrongly
  # typed matrix is rejected before the annealing, not hours into it.
  if (!is.integer(originalGenotypeMatrix)) {
    stop("originalGenotypeMatrix must already be stored as an integer matrix")
  }

  space <- prepareControlsSpace(genotypeMatrix, SVDReference, controlsMean,
                                casesMean, casesPDs, caseVariants,
                                controlsClustering)
  # caseVariants and caseCounts are consumed by different stages, so their
  # agreement is this function's to check -- again, before the annealing runs.
  if (nrow(as.matrix(caseCounts)) != length(space$variantRows)) {
    stop("caseCounts must have one row per case variant")
  }

  match <- function(candidates) {
    matchControlCandidates(originalGenotypeMatrix, caseCounts, candidates,
                           space$clusterIds, space$variantRows,
                           minLambda = minLambda,
                           softMinLambda = softMinLambda,
                           softMaxLambda = softMaxLambda,
                           maxLambda = maxLambda, min = min,
                           minCallRate = minCallRate)
  }

  # The sample names come from the matrix the counts were taken from.
  runControlSelection(space, match, colnames(originalGenotypeMatrix),
                      returnClusters, min, max, step, iterations, seed,
                      saThreads, exactPrecomputeThreads, exactClusterTileSize)
}

#' Select a set of controls without holding the control genotypes in memory.
#'
#' Does exactly what \code{\link{selectControls}} does, but reads the control
#' genotypes from files written by \code{\link{writeControlsMatrix}} instead of
#' taking them as matrices. The projection streams the reduced genotypes one
#' variant at a time and the matching streams the raw ones, so at no point does
#' a control genotype matrix exist in R. What is held is the reduced
#' \code{points} (one small vector per sample) and the per-cluster allele counts,
#' both of which are orders of magnitude smaller than the genotypes they come
#' from.
#'
#' \code{originalFile} may be a raw genotype matrix written with
#' \code{dtype = "u8"}, or the per-cluster counts produced by
#' \code{\link{collapseControlsMatrix}}. The collapsed form is both faster and
#' free of individual genotypes; because it has no samples left in it, the
#' selected controls come back as cluster identifiers.
#'
#' The two files are matched by name, not by position: the case variants are
#' resolved against the variant names of each file separately, so the files may
#' order their variants differently. The samples, however, must line up, and
#' with a collapsed file the clustering stored in it must be the one
#' \code{controlsClustering} describes. That last check compares the cluster
#' labels and the per-cluster sample counts, which is all a collapsed file
#' retains: a clustering that permutes samples between equally sized clusters of
#' the same name cannot be told apart from the one the file was built with.
#' @param genotypeFile path of the control matrix file holding the imputed
#' genotypes, written with \code{dtype} \code{"f64"} or \code{"f32"} (or
#' \code{"u8"} when the genotypes are integral).
#' @param originalFile path of the control matrix file holding the raw genotypes
#' (\code{dtype = "u8"}) or the collapsed per-cluster counts.
#' @param projectionThreads optional integer number of threads for the streaming
#' projection. By default uses all available hardware threads.
#' @inheritParams selectControls
#' @return the same list \code{\link{selectControls}} returns.
#' @seealso \code{\link{writeControlsMatrix}},
#' \code{\link{collapseControlsMatrix}}, \code{\link{selectControls}}
#' @export
selectControlsFromFiles <- function(genotypeFile, originalFile, casesPDs,
                                    casesMean, SVDReference, controlsMean,
                                    caseCounts, caseVariants = NULL,
                                    controlsClustering = NULL, minLambda = 0.75,
                                    softMinLambda = 0.9, softMaxLambda = 1.05,
                                    maxLambda = 1.3, min = 500, max = 1000,
                                    step = 50, iterations = 100000,
                                    minCallRate = 0.98, saThreads = NULL,
                                    exactPrecomputeThreads = NULL,
                                    exactClusterTileSize = 32L,
                                    projectionThreads = NULL,
                                    returnClusters = FALSE, seed = NULL) {
  seed <- resolveSeed(seed)
  returnClusters <- isTRUE(returnClusters)
  genotypeFile <- path.expand(genotypeFile)
  originalFile <- path.expand(originalFile)

  space <- prepareControlsSpace(genotypeFile, SVDReference, controlsMean,
                                casesMean, casesPDs, caseVariants,
                                controlsClustering, threads = projectionThreads)
  if (nrow(as.matrix(caseCounts)) != length(space$variantRows)) {
    stop("caseCounts must have one row per case variant")
  }

  info <- readControlsMatrixInfo(originalFile)
  if (info$dtype == "clusterCounts") {
    if (!identical(info$colnames, space$clusterLabels) ||
        !identical(as.integer(info$colweights),
                   as.integer(tabulate(space$clusterIds, length(info$colnames))))) {
      stop(paste("the cluster counts file was collapsed over a different",
                 "clustering than the one prepareControlsSpace resolved"))
    }
    # No samples survive in a collapsed file, so the selection can only be
    # reported as cluster identifiers.
    sampleNames <- NULL
  } else {
    if (!identical(info$colnames, space$sampleNames)) {
      stop(paste("genotypeFile and originalFile must hold the same samples in",
                 "the same order"))
    }
    sampleNames <- info$colnames
  }

  match <- function(candidates) {
    matchControlCandidates(originalFile, caseCounts, candidates,
                           clusterIds = space$clusterIds,
                           caseVariants = space$caseVariants,
                           minLambda = minLambda,
                           softMinLambda = softMinLambda,
                           softMaxLambda = softMaxLambda,
                           maxLambda = maxLambda, min = min,
                           minCallRate = minCallRate)
  }

  runControlSelection(space, match, sampleNames, returnClusters, min, max, step,
                      iterations, seed, saThreads, exactPrecomputeThreads,
                      exactClusterTileSize)
}
