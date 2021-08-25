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

matchControlsCluster <- function(cases, gmatrix, original, ...) {
  softMinLambda <- list(...)$softMinLambda
  softMaxLambda <- list(...)$softMaxLambda
  
  checkCaseInfo(cases)
  
  if (nrow(cases$counts) != nrow(gmatrix)){
    userError("Something is wrong with SNP selection")
  }
  
  results <- selectControls(gmatrix, original, cases$US, cases$mean, 
                            cases$counts, ...)
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
         clusters = setNames(good, cases$id), 
         ncontrols = if(good) length(results$controls) else 0)
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

countGoodControls <- function(l, r) {
  table <- rbind(l$table, r$table)
  length(unique(table$sample[table$cluster %in% goodClusters(l, r)]))
}

jointResult <- function(l, r, res) {
  ret <- list()
  table <- rbind(l$table, r$table)
  ret$table <- table
  ret$pvals <- c(l$pvals, r$pvals)
  ret$lambdas <- c(l$lambdas, r$lambdas)
  ret$ncontrols <- countGoodControls(l, r)
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

mergeCondition <- function(l, r, merged, mergeCoef) {
  if (subtreeFailed(l) && subtreeFailed(r)) {
    merged$minL < l$minL && merged$minL < r$minL
  } else {
    jointCount <- countGoodControls(l, r)
    mergeCoef * merged$ncontrols >  jointCount
  }
}

mergedOrJoint <- function(gmatrix, original, left, right, res, mergeCoef, 
                          ...) {
  if (mergeCondition(left, right, res, mergeCoef)) {  
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
    jointResult(left, right, res)
  }
}

recSelect <- function(gmatrix, original, cases, 
                      hierNode, clusterMergeCoef, ...) {
  cluster <- cases$population[[hierNode$id]]
  cluster$variants <- cases$variants
  res <- matchControlsCluster(cluster, gmatrix, original, ...)
  if (hierNode$type == "leaf") {
    return(res)
  } else {
    left <- recSelect(gmatrix, original, cases, hierNode$left,  
                      clusterMergeCoef, ...)
    right <- recSelect(gmatrix, original, cases, hierNode$right, 
                       clusterMergeCoef, ...)
    mergedOrJoint(gmatrix, original, left, right, res, clusterMergeCoef, ...)
  }
}

filter_variants <- function(population, ids) {
  population
}

#' Select a set of controls that matches a set of cases
#' @param controlGMatrix numeric matrix(0 - ref, 1 - het, 2 - both alt). 
#' Intermediate values are allowed, NAs are not.
#' @param cases result of calling function readInstanceFromYml.
#' @param originalControlGMatrix integer matrix(0 - ref, 1 - het, 2 - both alt).
#' Missing values are allowed.
#' @param clusterMergeCoef numeric coefficient of preference of merging clusters.
#' @param ... parameters to be passed to selectControls function.
#' @inheritParams selectControls
#' @export
selectControlsHier <- function(controlGMatrix, cases, 
                               originalControlGMatrix, 
                               clusterMergeCoef = 1.1, 
                               softMinLambda = 0.9, softMaxLambda = 1.05, 
                               ...) {
  stopifnot(all(!is.na(controlGMatrix)))
  stopifnot(all(rownames(controlGMatrix) == rownames(originalControlGMatrix)))
  
  #check if all sites from cases are found in controls
  
  extra_variants <- which(!(cases$variants %in% rownames(controlGMatrix)))
  if (length(extra_variants) != 0) {
    num <- min(5, length(extra_variants))
    userError(paste("Not all case sites are found in controls. Please use",
               "list of available sites to create your case genotype",
               "matrix. These variants must not be present: ", 
               paste(cases$variants[extra_variants][1:num], 
                     collapse = ", ")))
  }
  
  variants <- intersect(cases$variants, rownames(controlGMatrix))
  controlGMatrix <- controlGMatrix[variants, ]
  originalControlGMatrix <- originalControlGMatrix[variants, ]
  
  cases_variants_ids <- setNames(1:length(variants), variants)[cases$variants] 
  cases$population <- lapply(cases$population, function(population) {
    population$counts <- population$counts[cases_variants_ids, ]
    population
  })
  
  recSelect(controlGMatrix, originalControlGMatrix,  
            cases, cases$hierarchy, clusterMergeCoef, 
            softMinLambda = softMinLambda, softMaxLambda = softMaxLambda, ...)
}
