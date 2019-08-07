condition <- function(subclass, message, call = sys.call(-1), ...) {
  structure(
    class = c(subclass, "condition"),
    list(message = message, call = call),
    ...
  )
}

checkCaseInfo <- function(cases, variants) {
  if (nrow(cases$counts) != nrow(cases$US)) {
    userError("Case counts table does not correspond to U matrix")
  }
  if (!is.numeric(cases$US) | any(is.na(cases$US))) {
    userError("Matrix U contains missing data or non-numeric values")
  }
  
  ##check if all sites from cases are found in controls
  
  extra_variants <- which(!(cases$variants %in% variants))
  if (length(extra_variants) != 0) {
    num <- min(5, length(extra_variants))
    userError(paste("Not all case sites are found in controls. Please use",
               "list of available sites to create your case genotype",
               "matrix. These variants must not be present: ", 
               paste(cases$variants[extra_variants][1:num], 
                     collapse = ", ")))
  }
}

matchControlsCluster <- function(cases, gmatrix, original, variants, ...) {
  softMinLambda <- list(...)$softMinLambda
  softMaxLambda <- list(...)$softMaxLambda
  print(cases$id)
  
  checkCaseInfo(cases, variants)
  sharedSites <- intersect(cases$variants, variants)
  selector <- which(variants %in% sharedSites)
  gmatrix <- gmatrix[selector, ]
  original <- original[selector, ]
  
  if (nrow(cases$counts) != nrow(gmatrix)){
    userError("Something is wrong with SNP selection")
  }
  
  U <- apply(cases$US, 2, function(x) x / sqrt(sum(x ^ 2))) 
  results <- selectControls(gmatrix, original, U, cases$counts, ...)
  resid <- results$residuals[results$controls]
  df <- data.frame(sample = names(resid), value = resid, 
                   cluster = if (length(resid) == 0) c() else cases$id,
                   stringsAsFactors = FALSE, row.names = NULL)
  pvals <- list()
  pvals[[as.character(cases$id)]] <- results$pvals
  if (length(results$controls) > 0) {
    lam <- results$optimal_lambda
    good <- lam > softMinLambda && lam < softMaxLambda
    list(table = df, pvals = pvals, minL = results$optimal_lambda, cases = cases, 
         clusters = setNames(good, cases$id), 
         ncontrols = if(good) length(results$controls) else 0)
  } else {
    list(table = data.frame(), pvals = c(), minL = Inf, cases = cases,
         clusters = c(), ncontrols = 0)
  }
}

distributeControls <- function(residuals, threshold = 100) {
  if (nrow(residuals) == 0) {
    return(NULL)
  }
  residuals <- residuals[order(residuals$value), ]
  while (TRUE) {
    if (nrow(residuals) == 0) {
      return(residuals)
    }
    top <- residuals[!duplicated(residuals$sample), ]
    clusterCounts <- table(top$cluster)
    outsider <- which.min(clusterCounts)
    if (clusterCounts[outsider] >= threshold) {
      return(top)
    } else {
      residuals <- residuals[residuals$cluster != outsider, ]
    }
  }
  residuals
}

mergeCases <- function(left, right) {
  ret <- list()
  ret$id <- paste0(left$id, ", ", right$id)
  ret$counts <- left$counts + right$counts
  lUS <- left$US
  rUS <- right$US
  svd <- RSpectra::svds(cbind(lUS, rUS), k = max(ncol(lUS), ncol(rUS)))
  ret$US <- svd$u %*% diag(svd$d)
  ret$variants <- left$variants
  ret
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

jointResult <- function(l, r) {
  ret <- list()
  table <- rbind(l$table, r$table)
  ret$table <- table
  ret$pvals <- c(l$pvals, r$pvals)
  goodClusters <- goodClusters(l, r)
  ret$ncontrols <- countGoodControls(l, r)
  ret$clusters <- c(l$clusters, r$clusters)
  ret$minL <- min(l$minL, r$minL)
  if (l$ncontrols == 0 && r$ncontrols > 0) {
    ret$cases <- r$cases
  } else if (r$ncontrols == 0 && l$ncontrols > 0) {
    ret$cases <- l$cases
  } else {
    ret$cases <- mergeCases(l$cases, r$cases)
  }
  ret
}

subtreesFailed <- function(l, r) {
  return(l$ncontrols == 0 && r$ncontrols == 0)
}

mergeCondition <- function(l, r, merged, mergeCoef) {
  if (subtreesFailed(l, r)) {
    merged$minL < l$minL && merged$minL < r$minL
  } else {
    jointCount <- countGoodControls(l, r)
    mergeCoef * merged$ncontrols >  jointCount
  }
}

mergedOrJoint <- function(gmatrix, original, variants, left, right, mergeCoef, 
                          ...) {
  cases <- mergeCases(left$cases, right$cases)
  res <- matchControlsCluster(cases, gmatrix, original, variants, ...)
  if (mergeCondition(left, right, res, mergeCoef)) {  
    if (!subtreesFailed(left, right)) {
      cls <- goodClusters(left, right)
      table <- rbind(left$table, right$table)
      table <- table[!(table$cluster %in% cls), ]
    } else {
      table <- NULL
    }
    res$table <- rbind(table, res$table)
    res
  } else {
    jointResult(left, right)
  }
}

recSelect <- function(gmatrix, original, variants, cases, 
                      hierNode, clusterMergeCoef, ...) {
  if (hierNode$type == "leaf") {
    cluster <- cases$population[[hierNode$id]]
    cluster$variants <- cases$variants
    matchControlsCluster(cluster, gmatrix, original, variants, ...)
  } else {
    left <- recSelect(gmatrix, original, variants, cases, hierNode$left,  
                      clusterMergeCoef, ...)
    right <- recSelect(gmatrix, original, variants, cases, hierNode$right, 
                       clusterMergeCoef, ...)
    mergedOrJoint(gmatrix, original, variants, left, right, clusterMergeCoef, ...)
  }
}

#' Select a set of controls that matches a set of cases
#' @param controlGMatrix numeric matrix(0 - ref, 1 - het, 2 - both alt). 
#' Intermediate values are allowed, NAs are not.
#' @param controlVariants character vector of variants(format: "chrN:POS\\tREF\\tALT") 
#' @param cases result of calling function readInstanceFromYml.
#' @param imputationMatrix logical matrix of the same dimentions as 
#' controlGMatrix. Element is set to TRUE if the corresponding value has been
#' imputed in controlGMatrix.
#' @param originalControlGMatrix integer matrix(0 - ref, 1 - het, 2 - both alt).
#' Missing values are allowed.
#' @param clusterMergeCoef numeric coefficient of preference of merging clusters.
#' @param ... parameters to be passed to selectControls function.
#' @inheritParams selectControls
#' @export
selectControlsHier <- function(controlGMatrix, controlVariants, cases, 
                               imputationMatrix = NULL, 
                               originalControlGMatrix = NULL, 
                               clusterMergeCoef = 1.1, 
                               softMinLambda = 0.9, softMaxLambda = 1.05, 
                               ...) {
  stopifnot(all(!is.na(controlGMatrix)))
  if (is.null(originalControlGMatrix)) {
    originalControlGMatrix <- controlGMatrix
    mode(originalControlGMatrix) <- "integer"
    if (!is.null(imputationMatrix)) {
      originalControlGMatrix[imputationMatrix] <- NA
    }
  }
  
  recSelect(controlGMatrix, originalControlGMatrix, controlVariants, 
            cases, cases$hierarchy, clusterMergeCoef, 
            softMinLambda = softMinLambda, softMaxLambda = softMaxLambda, ...)
}