checkCaseInfo <- function(cases, variants, min_controls) {
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
  
  if (length(cases$variants) < min_controls) {
    stop("Too little variants available for analysis.")
  }
}

matchControlsCluster <- function(cases, gmatrix, original, variants,
                                 controlsThreshold, ...) {
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
    resid <- results$residuals
    data.frame(sample = names(resid), value = resid, cluster = cases$cluster, 
               stringsAsFactors = FALSE, row.names = NULL)
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
  ret$cluster <- paste(left$cluster, " & ", right$cluster)
  ret$counts <- left$counts + right$counts
  svd <- RSpectra::svds(cbind(left$US, right$US), k = max(ncol(left$US), 
                                                          ncol(right$US)))
  ret$US <- svd$u %*% diag(svd$d)
  ret
}

recSelect <- function(gmatrix, original, variants, cases, 
                      hierNode, threshold, separate, ...) {
  if (hierNode$type == "leaf") {
    curr <- cases[[hierNode$id]]
    curr$variants <- cases$variants
    controls <- matchControlsCluster(curr, gmatrix, original, variants, ...)
    list(cases = curr, controls = unique(controls$sample), table = controls)
  } else {
    left <- recSelect(gmatrix, original, variants, cases, hierNode$left,  ...)
    right <- recSelect(gmatrix, original, variants, cases, hierNode$right, ...)
    mergedCases <- mergeCases(left$cases, right$cases)
    merged <- matchControlsCluster(mergeCases(left$cases, righ$cases), gmatrix, 
                                  original, variants, threshold, ...)
    jointControls <- union(left$controls, right$controls)
    if (length(merged$controls) >= length(jointControls)) {
      list(cases = mergedCases, controls = merged$controls, table = merged)
    } else {
      table <- rbind(left$table, right$table)
      if (separate) {
        table <- distributeControls(table)
        jointControls <- unique(table$sample)
      }
      list(cases = mergedCases, controls = jointControls, table = table)
    }
  }
}

#' Select a set of controls that matches a set of cases
#' @param controlGMatrix numeric matrix(0 - ref, 1 - het, 2 - both alt). 
#' Intermediate values are allowed, NAs are not.
#' @param controlVariants character vector of variants(format: "chrN:POS\\tREF\\tALT") 
#' @param  cases result of calling function readInstanceFromYml.
#' @param  imputationMatrix logical matrix of the same dimentions as 
#' controlGMatrix. Element is set to TRUE if the corresponding value has been
#' imputed in controlGMatrix.
#' @param originalControlGMatrix integer matrix(0 - ref, 1 - het, 2 - both alt).
#' Missing values are allowed.
#' @param separate logical if TRUE results will not contain overlapping 
#' sets of controls.
#' @param ... parameters to be passed to selectControls function.
#' @export
selectControlsHier <- function(controlGMatrix, controlVariants, cases, 
                               imputationMatrix = NULL, 
                               originalControlGMatrix = NULL, 
                               separate = FALSE, minControls, ...) {
  stopifnot(all(!is.na(controlGMatrix)))
  if (is.null(originalControlGMatrix)) {
    originalControlGMatrix <- controlGMatrix
    mode(originalControlGMatrix) <- "integer"
    if (!is.null(imputationMatrix)) {
      originalControlGMatrix[imputationMatrix] <- NA
    }
  }
  
  checkCaseInfo(cases$variants, variants)
  selection <- recSelect(controlGMatrix, originalControlGMatrix, 
                         cases$population, cases$hierarchy, minControls, 
                         separate, ...)
}