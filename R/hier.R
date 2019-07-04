checkCaseInfo <- function(cases, variants) {
  if (nrow(cases$counts) != nrow(cases$US)) {
    stop("Case counts table does not correspond to U matrix")
  }
  if (!is.numeric(cases$US) | any(is.na(cases$US))) {
    stop("Matrix U contains missing data or non-numeric values")
  }
  
  ##check if all sites from cases are found in controls
  
  extra_variants <- which(!(cases$variants %in% variants))
  if (length(extra_variants) != 0) {
    num <- min(5, length(extra_variants))
    stop(paste("Not all case sites are found in controls. Please use",
               "list of available sites to create your case genotype",
               "matrix. These variants must not be present: ", 
               paste(cases$variants[extra_variants][1:num], 
                     collapse = ", ")))
  }
}

matchControlsCluster <- function(cases, gmatrix, original, variants,
                                 controlsThreshold, softMinLambda, 
                                 softMaxLambda, ...) {
  checkCaseInfo(cases, variants)
  sharedSites <- intersect(cases$variants, variants)
  selector <- which(variants %in% sharedSites)
  gmatrix <- gmatrix[selector, ]
  original <- original[selector, ]
  
  if (nrow(cases$counts) != nrow(gmatrix)){
    stop("Something is wrong with SNP selection")
  }
  
  U <- apply(cases$US, 2, function(x) x / sqrt(sum(x ^ 2))) 
  results <- selectControls(gmatrix, original, U, cases$counts, ...)
  
  resid <- results$residuals[results$controls]
  df <- data.frame(sample = names(resid), value = resid, 
                   cluster = cases$cluster, stringsAsFactors = FALSE, 
                   row.names = NULL)
  pvals <- list()
  pvals[[cases$cluster]] <- results$pvals
  if (length(results$controls) >= controlsThreshold) {
    lam <- results$optimalLambda
    good <- lam > softMinLambda && lam < softMaxLambda
    list(table = df, pvals = pvals, minL = results$optimalLambda, cases = cases, 
         clusters = setNames(good, cases$cluster), 
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
  ret$cluster <- paste0(left$cluster, ", ", right$cluster)
  ret$counts <- left$counts + right$counts
  svd <- RSpectra::svds(cbind(left$US, right$US), k = max(ncol(left$US), 
                                                          ncol(right$US)))
  ret$US <- svd$u %*% diag(svd$d)
  ret$variants <- left$variants
  ret
}

goodClusters <- function(l, r) {
  c(l$clusters[l$clusters], r$clusters[r$clusters])
}

countGoodControls <- function(l, r) {
  length(unique(table$sample[table$cluster %in% goodClusters]))
}

jointResult <- function(l, r) {
  ret <- list()
  table <- rbind(l$table, r$table)
  ret$table <- table
  ret$pvals <- c(l$pvals, r$pvals)
  goodClusters <- goodClusters(l, r)
  ret$ncontrols <- countGoodControls(l, r)
  ret$clusters <- c(left$clusters, right$clusters)
  if (l$ncontrols == 0 && r$ncontrols > 0) {
    ret$cases <- r$cases
  } else if (r$ncontrols == 0 && l$ncontrols > 0) {
    ret$cases <- l$cases
  } else {
    ret$cases <- mergeCases(l$cases, r$cases)
  }
  ret
}

mergeCondition <- function(l, r, merged, mergeCoef) {
  if (l$ncontrols == 0 && r$ncontrols == 0) {
    merged$minL < l$minL && merged$minL < r$minL
  } else {
    jointCount <- countGoodControls(l, r)
    mergeCoef * merged$ncontrols >  jointCount
  }
}

mergedOrJoint <- function(gmatrix, original, variants, left, right, mergeCoef, 
                          ...) {
  cases <- mergeCases(left, right)
  res <- matchControlsCluster(cases$cluster, gmatrix, original, variants, ...)
  if (mergeCondition(l, r, res, mergeCoef)) {  
    cls <- goodClusters(l, r)
    table <- rbind(l$table, r$table)
    table <- table[, !(table$cluster %in% cls)]
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
                      clusterMergeCoef)
    right <- recSelect(gmatrix, original, variants, cases, hierNode$right, 
                       clusterMergeCoef)
    mergedOrJoint(left, right)
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
#' @param minControls integer minimum number of controls selected for 
#' each cluster
#' @param clusterMergeCoef numeric coefficient of preference of merging clusters.
#' @param ... parameters to be passed to selectControls function.
#' @inheritParams selectControls
#' @export
selectControlsHier <- function(controlGMatrix, controlVariants, cases, 
                               imputationMatrix = NULL, 
                               originalControlGMatrix = NULL, 
                               minControls = 50, clusterMergeCoef = 1.1, 
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
            cases, cases$hierarchy, minControls, clusterMergeCoef, 
            softMinLambda, softMaxLambda, ...)
}