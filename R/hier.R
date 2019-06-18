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
                                 controlsThreshold, ...) {
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
  
  if (length(results$controls) > controlsThreshold) {
    resid <- results$residuals[results$controls]
    df <- data.frame(sample = names(resid), value = resid, 
                     cluster = cases$cluster, stringsAsFactors = FALSE, 
                     row.names = NULL)
    list(resid = df, pvals = results$pvals, cl = cases$cluster, 
         lambda = results$optimal_lambda)
  } else {
    NULL
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

recSelect <- function(gmatrix, original, variants, cases, 
                      hierNode, threshold, clusterMergeCoef, 
                      softMinLambda, softMaxLambda, 
                      ...) {
  filterSamples <- function(l, samples) {
    if(l >= softMinLambda && l <= softMaxLambda) {
      samples
    } else {
      c()
    }
  }
  if (hierNode$type == "leaf") {
    curr <- cases$population[[hierNode$id]]
    curr$variants <- cases$variants
    matching <- matchControlsCluster(curr, gmatrix, original, variants, 
                                     threshold, ...)
    controls <- matching$resid
    pvals <- list()
    pvals[[matching$cl]] <- matching$pvals
    samples <- filterSamples(matching$lambda, controls$sample)
    
    list(cases = curr, controls = samples, table = controls, pvals = pvals)
  } else {
    left <- recSelect(gmatrix, original, variants, cases, hierNode$left,  
                      threshold, clusterMergeCoef,  softMinLambda,
                      softMaxLambda, ...)
    right <- recSelect(gmatrix, original, variants, cases, hierNode$right, 
                       threshold, clusterMergeCoef, softMinLambda, 
                       softMaxLambda, ...)
    mergedCases <- mergeCases(left$cases, right$cases)
    mergedMatching <- matchControlsCluster(mergeCases(left$cases, right$cases), 
                                           gmatrix, original, variants, 
                                           threshold, ...)
    mergedSample <- filterSamples(mergedMatching$lambda, mergedMatching$resid$sample)
    merged <- mergedMatching$resid
    jointControls <- union(left$controls, right$controls)
    if (length(mergedSample) >= (1.0 / clusterMergeCoef) * length(jointControls)) {
      pvals <- list()
      pvals[[mergedMatching$cl]] <- mergedMatching$pvals
      list(cases = mergedCases, controls = merged$sample, table = merged, 
           pvals = pvals)
    } else {
      table <- rbind(left$table, right$table)
      pvals <- c(left$pvals, right$pvals)
      list(cases = mergedCases, controls = jointControls, table = table, 
           pvals = pvals)
    }
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