condition <- function(subclass, message, call = sys.call(-1), ...) {
  structure(
    class = c(subclass, "condition"),
    list(message = message, call = call),
    ...
  )
}

checkCaseInfo <- function(cases) {
  if (!is.numeric(cases$US) | any(is.na(cases$US))) {
    userError("Matrix U contains missing data or non-numeric values")
  }
}

# What one returned control entry is worth: one sample, or the size of the
# cluster when the selection comes back as cluster identifiers. Without this the
# merge heuristic would weigh cluster counts against sample counts.
controlsWeight <- function(controls, controlWeights) {
  if (is.null(controlWeights)) {
    return(length(controls))
  }
  sum(controlWeights[controls])
}

matchControlsCluster <- function(cases, selector, controlWeights, ...) {
  softMinLambda <- list(...)$softMinLambda
  softMaxLambda <- list(...)$softMaxLambda
  
  checkCaseInfo(cases)
  
  if (nrow(cases$counts) != length(cases$variants)){
    userError("Something is wrong with SNP selection")
  }
  
  results <- selector(cases)
  df <- data.frame(sample = results$controls, 
                   cluster = if (length(results$controls) == 0) c() else cases$id,
                   stringsAsFactors = FALSE, row.names = NULL)
  pvals <- list()
  lambdas <- list()
  clust <- as.character(cases$id) 
  pvals[[clust]] <- results$pvals
  lambdas[[clust]] <- results$optimal_lambda
  if (length(results$controls) > 0) {
    lam <- results$optimal_lambda
    good <- lam > softMinLambda && lam < softMaxLambda
    list(table = df, pvals = pvals, lambdas = lambdas, 
         minL = results$optimal_lambda, cases = cases, 
         clusters = stats::setNames(good, cases$id), 
         ncontrols = if(good) controlsWeight(results$controls, controlWeights) else 0)
  } else {
    list(table = data.frame(), pvals = c(), lambdas = c(), minL = Inf, 
         cases = cases, clusters = c(), ncontrols = 0)
  }
}

goodClusters <- function(l, r) {
  ret <- names(c(l$clusters[l$clusters], r$clusters[r$clusters]))
  if (is.null(ret)) {
    ret <- character(0)
  }
  ret
}

countGoodControls <- function(l, r, controlWeights) {
  table <- rbind(l$table, r$table)
  controlsWeight(unique(table$sample[table$cluster %in% goodClusters(l, r)]),
                 controlWeights)
}

jointResult <- function(l, r, res, controlWeights) {
  ret <- list()
  table <- rbind(l$table, r$table)
  ret$table <- table
  ret$pvals <- c(l$pvals, r$pvals)
  ret$lambdas <- c(l$lambdas, r$lambdas)
  ret$ncontrols <- countGoodControls(l, r, controlWeights)
  ret$clusters <- c(l$clusters, r$clusters)
  ret$minL <- min(l$minL, r$minL)
  if (l$ncontrols == 0 && r$ncontrols > 0) {
    ret$cases <- r$cases
  } else if (r$ncontrols == 0 && l$ncontrols > 0) {
    ret$cases <- l$cases
  } else {
    ret$cases <- res$cases
  }
  ret
}

subtreeFailed <- function(t) {
  t$ncontrols == 0
}

mergeCondition <- function(l, r, merged, mergeCoef, controlWeights) {
  if (subtreeFailed(l) && subtreeFailed(r)) {
    merged$minL < l$minL && merged$minL < r$minL
  } else {
    jointCount <- countGoodControls(l, r, controlWeights)
    mergeCoef * merged$ncontrols >  jointCount
  }
}

mergedOrJoint <- function(left, right, res, mergeCoef, controlWeights, ...) {
  if (mergeCondition(left, right, res, mergeCoef, controlWeights)) {  
    cls <- goodClusters(left, right)
    if (subtreeFailed(left)) {
      cls <- c(cls, names(left$clusters))
    }
    if (subtreeFailed(right)) {
      cls <- c(cls, names(right$clusters))
    }
    table <- rbind(left$table, right$table)
    table <- table[!(table$cluster %in% cls), ]
    
    allClusters <- c(left$clusters, right$clusters, res$clusters)
    res$pvals <- c(left$pvals, right$pvals, res$pvals)
    res$lambdas <- c(left$lambdas, right$lambdas, res$lambdas)
    keep <- setdiff(names(allClusters), cls)
    res$clusters <- allClusters[keep]
    for (cl in names(res$pvals)) {
      if (!(cl %in% keep)) {
        res$pvals[[cl]] <- NULL
        res$lambdas[[cl]] <- NULL
      }
    }
    res$table <- rbind(table, res$table)
    res
  } else {
    jointResult(left, right, res, controlWeights)
  }
}

recSelect <- function(selector, cases, hierNode, clusterMergeCoef, controlWeights,
                      ...) {
  cluster <- cases$population[[hierNode$id]]
  cluster$variants <- cases$variants
  res <- matchControlsCluster(cluster, selector, controlWeights, ...)
  if (hierNode$type == "leaf") {
    return(res)
  } else {
    left <- recSelect(selector, cases, hierNode$left, clusterMergeCoef,
                      controlWeights, ...)
    right <- recSelect(selector, cases, hierNode$right, clusterMergeCoef,
                       controlWeights, ...)
    mergedOrJoint(left, right, res, clusterMergeCoef, controlWeights, ...)
  }
}

filter_variants <- function(population, ids) {
  population
}

# Everything the hierarchy needs that does not depend on which case cluster is
# being matched: the reduced control coordinates and the per-cluster counts. Both
# are built once and re-aimed per node, which is what keeps the walk over a
# hierarchy of H nodes from repeating H projections and H collapses.
hierSelector <- function(space, counts, sampleNames, dots) {
  function(cluster) {
    do.call(selectWithPrepared,
            c(list(space = space, counts = counts, sampleNames = sampleNames,
                   casesMean = cluster$mean, casesPDs = cluster$US,
                   caseCounts = cluster$counts), dots))
  }
}

hierControlWeights <- function(space, sampleNames, returnClusters) {
  if (!isTRUE(returnClusters) && !is.null(sampleNames)) {
    return(NULL)
  }
  stats::setNames(tabulate(space$clusterIds, length(space$clusterLabels)),
                  space$clusterLabels)
}

runHierSelection <- function(space, counts, sampleNames, cases, clusterMergeCoef,
                             dots) {
  selector <- hierSelector(space, counts, sampleNames, dots)
  controlWeights <- hierControlWeights(space, sampleNames, dots$returnClusters)
  ret <- do.call(recSelect,
                 c(list(selector, cases, cases$hierarchy, clusterMergeCoef,
                        controlWeights), dots))
  ret$names <- cases$names
  ret
}

#' Select a set of controls that populationally matches a set of cases.
#'
#' Walks the case hierarchy and matches every node against the same controls.
#' The control side of the work -- projecting the controls into the case
#' PCA-like space and collapsing them into per-cluster allele counts -- does not
#' depend on which case cluster is being matched, so it is done once and reused
#' for every node; only the multivariate normality search and the lambda scoring
#' run per node.
#' @param controlGMatrix numeric matrix(0 - ref, 1 - het, 2 - both alt).
#' Intermediate values are allowed, NAs are not. Rows are named by variant and
#' must include every entry of \code{cases$variants}; extra control variants are
#' allowed and are matched by name through an index view rather than by slicing.
#' @param originalControlGMatrix integer matrix(0 - ref, 1 - het, 2 - both alt)
#' with missing values allowed. Its rows must correspond to those of
#' \code{controlGMatrix}.
#' @param cases result of calling function readInstanceFromYml.
#' @param clusterMergeCoef numeric coefficient of preference of merging clusters.
#' @param ... parameters to be passed to selectControls function.
#' @inheritParams selectControls
#' @seealso \code{\link{selectControlsHierFromFiles}}
#' @export
selectControlsHier <- function(controlGMatrix, originalControlGMatrix, 
                               cases, SVDReference, controlsMean, 
                               clusterMergeCoef = 1.1, 
                               softMinLambda = 0.9, softMaxLambda = 1.05, 
                               ...) {
  stopifnot(all(!is.na(controlGMatrix)))
  stopifnot(all(rownames(controlGMatrix) == rownames(originalControlGMatrix)))

  if (!all(cases$variants %in% rownames(controlGMatrix))) {
    userError("cases$variants must be a subset of the control genotype matrix rows")
  }

  dots <- c(list(softMinLambda = softMinLambda, softMaxLambda = softMaxLambda),
            list(...))
  space <- hierControlsSpace(controlGMatrix, cases, SVDReference, controlsMean, dots)
  counts <- prepareClusterCounts(originalControlGMatrix, space$clusterIds,
                                 space$variantRows)
  runHierSelection(space, counts, colnames(originalControlGMatrix), cases,
                   clusterMergeCoef, dots)
}

# Any population will do to seed the space: casesMean and casesPDs only set the
# distribution the search aims at, and every node re-aims it.
hierControlsSpace <- function(controls, cases, SVDReference, controlsMean, dots) {
  seed <- cases$population[[1]]
  prepareControlsSpace(controls, SVDReference, controlsMean, seed$mean, seed$US,
                       caseVariants = cases$variants,
                       controlsClustering = dots$controlsClustering,
                       threads = dots$projectionThreads)
}

#' Match a case hierarchy without holding the control genotypes in memory.
#'
#' Does what \code{\link{selectControlsHier}} does, reading the controls from
#' files written by \code{\link{writeControlsMatrix}} instead of taking them as
#' matrices, exactly as \code{\link{selectControlsFromFiles}} does for a single
#' set of cases. Because the control side is prepared once for the whole
#' hierarchy, each file is read once no matter how many nodes the hierarchy has.
#'
#' With a collapsed counts file the selected controls come back as cluster
#' identifiers, and the merge heuristic weighs each of them by the number of
#' samples its cluster holds.
#' @param genotypeFile path of the control matrix file holding the imputed
#' genotypes.
#' @param originalFile path of the control matrix file holding the raw genotypes
#' (\code{dtype = "u8"}) or the collapsed per-cluster counts.
#' @param projectionThreads optional integer number of threads for the streaming
#' projection. By default uses all available hardware threads.
#' @inheritParams selectControlsHier
#' @inheritParams selectControls
#' @seealso \code{\link{selectControlsHier}}, \code{\link{writeControlsMatrix}},
#' \code{\link{collapseControlsMatrix}}
#' @export
selectControlsHierFromFiles <- function(genotypeFile, originalFile, cases,
                                        SVDReference, controlsMean,
                                        clusterMergeCoef = 1.1,
                                        softMinLambda = 0.9, softMaxLambda = 1.05,
                                        projectionThreads = NULL, ...) {
  genotypeFile <- path.expand(genotypeFile)
  originalFile <- path.expand(originalFile)

  dots <- c(list(softMinLambda = softMinLambda, softMaxLambda = softMaxLambda,
                 projectionThreads = projectionThreads), list(...))
  space <- hierControlsSpace(genotypeFile, cases, SVDReference, controlsMean, dots)
  sampleNames <- controlsFileSampleNames(originalFile, space)
  counts <- prepareClusterCounts(originalFile, clusterIds = space$clusterIds,
                                 caseVariants = space$caseVariants)
  runHierSelection(space, counts, sampleNames, cases, clusterMergeCoef, dots)
}
