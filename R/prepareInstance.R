#' Estimate dendrogram for case clustering
#' uses combiM output of clustCombi() function
#' from Mclust package
#' @param combiM combiM output of clustCombi() function
dendrogramEstimate <- function(combiM){
  stopifnot(class(combiM) == "list")
  G <- length(combiM)
  curr <- 1:G
  merged <- -(1:G)
  merge <- matrix(NA, G - 1, 2)
  for (k in 1:(G - 1)) {
    Kp = G - k + 1
    l1 = which(!combiM[[Kp - 1]] %*% rep(1, Kp) == 1)
    l2 = (combiM[[Kp - 1]] %*% curr)[l1] - curr[l1]
    curr <- setdiff(curr, max(l1, l2))
    merge[k, ] <- merged[c(l1, l2)]
    merged[merged == merged[l1] | merged == merged[l2]] <- k
  }
  merge
}

collapsingToTree <- function(collapsing) {
  nodes <- list()
  if (is.null(collapsing)) {
    return(1)
  }
  
  numToNode <- function(x) {
    if(x < 0) {
      -x
    } else {
      nodes[[x]]
    }
  }
  
  for(i in 1:nrow(collapsing)) {
    nodes[[i]] <- list(numToNode(collapsing[i, 1]), numToNode(collapsing[i, 2]))
  } 
  nodes[[length(nodes)]]
}

#' Estimate clusters within case data.
#' 
#' Clusters are estimated from PC loadings using
#' Gaussian mixed model fitting from Mclust package
#' using ellipsoidal, varying volume, shape, and orientation
#' model ("VVV").
#' @param PCA results of gmatrix principal component analysis.
#' @param plotBIC logical whether a plot of Bayesian Information 
#' Criterion vs cluster number should be returned.
#' @param plotDendrogram whether dendrogram should be plotted.
#' @param minClusters minimal number of clusters.
#' @param clusters maximum number of clusters.
#' @param keepSamples names of the samples to keep in results. The feature is
#' used in order to cluster jointly with another dataset.
#' @export
estimateCaseClusters <- function (PCA, plotBIC = FALSE, plotDendrogram = FALSE, 
                                   minClusters = 1, clusters = 20, keepSamples = NULL) {
  stopifnot(class(PCA) %in% c("matrix", "data.frame", "array"))
  clResults <- mclust::Mclust(data = PCA, G = minClusters:clusters, modelNames = "VVV")
  if (is.null(clResults)) {
    stop("Couldn't fit gaussian mixed model to data. Try another range of clusters.") 
  }
  if (plotBIC) {
    mclust::plot.Mclust(x = clResults, what = "BIC", xlab = "Number of Components", 
                        ylab = "Bayesian Information Criterion")
  }
  clusters <- clResults$G
  if (clusters > 1) {
    clCombi <- mclust::clustCombi(clResults)
    collapsing <- dendrogramEstimate(clCombi$combiM)
  } else {
    collapsing <- NULL
  }
  if (plotDendrogram & clusters > 1) {
    mclust::combiTree(clCombi, type = "rectangle", yaxis = "entropy")
  }
  if(is.null(keepSamples)){
    keepSamples <- rownames(PCA)
  }
  res <- clustering(clResults$classification, collapsingToTree(collapsing))
  names(res$samples) <- rownames(PCA)
  if(!all(rownames(PCA) %in% keepSamples)){
    clustersToRemove <- setdiff(clResults$classification, clResults$classification[keepSamples])
    clustersToRemove <- unique(clustersToRemove)
    if(length(clustersToRemove) > 0){
      for(i in 1:length(clustersToRemove)){
        res <- SVDFunctions::removeCluster(res, clustersToRemove[i])
      }
    }
  }
  res$samples <- res$samples[keepSamples]
  res
}

#' Estimate PCA from gmatrix without missing values.
#' 
#' PCA is estimated using SVD of the mean centered genotype matrix.
#' PCA coordinates will be needed to perform clustering of the 
#' case samples to ensure relatively homogenious matching.
#' @param gmatrix Genotype matrix wihout missing values.
#' @param SVDReference a U matrix of SVD decomposition of centered reference
#' genotype matrix. 
#' @param referenceMean numeric vector representing per-variant genotype mean 
#' value of the genotype matrix that were used to construct \code{SVDReference}.
#' @param components Number of principal components to be computed.
#' @export
gmatrixPCA <- function(gmatrix, SVDReference = NULL, referenceMean = NULL, components = 10){
  if (is.null(referenceMean)) {
    referenceMean <- rowMeans(gmatrix)
  }
  
  if (is.null(SVDReference)) {
    SVDReference <- RSpectra::svds(gmatrix - referenceMean, k = components)$u
    rownames(SVDReference) <- rownames(gmatrix)
  }
  
  variants <- intersect(rownames(SVDReference), rownames(gmatrix))
  gmatrix <- gmatrix[variants, ]
  
  refVariants <- stats::setNames(1:nrow(SVDReference), rownames(SVDReference))
  ids <- refVariants[variants]
  
  SVDReference <- SVDReference[ids, ]
  referenceMean <- referenceMean[ids]
  
  gmatrix <- gmatrix - referenceMean
  pca <- t(SVDReference) %*% gmatrix
  
  pca <- pca[1:components, ]
  pca <- t(pca)
  colnames(pca) <- c(paste("PC", 1:components, sep = ""))
  rownames(pca) <- colnames(gmatrix)
  pca
}

#' Filter gmatrix by minimal sample call rate and minimal variant call rate.
#' 
#' @param minVariantCallRate minimal variant call rate.
#' @param minSampleCallRate minimal sample call rate.
#' @inheritParams prepareInstance
#' @export
filterGmatrix <- function(gmatrix, imputationResults, minVariantCallRate = 0.95,
                          minSampleCallRate = 0.95) {
  nSamples <- ncol(gmatrix)
  nVariants <- nrow(gmatrix)
  sampleCallRate <- apply(imputationResults, 2, function(x){
    sum(!x) / length(x)
  })
  gmatrix <- gmatrix[, which(sampleCallRate > minSampleCallRate)]
  
  variantCallRate <- apply(imputationResults, 1, function(x) {
    sum(!x) / length(x) 
  })
  
  gmatrix <- gmatrix[which(variantCallRate > minVariantCallRate), ]
  imputationResults <- imputationResults[variantCallRate > minVariantCallRate, 
                                         sampleCallRate > minSampleCallRate]
  
  message("Kept ", ncol(gmatrix), " out of ", nSamples, " individuals")
  message("Kept ", nrow(gmatrix), " out of ", nVariants, " variants")
  list(gmatrix = gmatrix, imputationResults = imputationResults)
}

#' Plot 3d PCA plot.
#' 
#' @param PCA pca of genotype matrix returned by function gmatrixPCA.
#' @inheritParams prepareInstance
#' @inheritParams estimateCaseClusters
#' @export 
plotPCA <- function(PCA, clusters = NULL) {
  if (!is.null(clusters)) {
    stopifnot(all(rownames(PCA) %in% names(clusters$samples)))
  }
  samples <- rownames(PCA)
  colors <- if (is.null(clusters)) 1 else clusters$classes[clusters$samples[samples]]
  texts <- paste("Sample:", rownames(PCA))
  if (!is.null(clusters)) {
    texts <- paste0(texts, "\nCluster: ", clusters$classes[clusters$samples[samples]])
  }
  plotly::add_markers(plotly::plot_ly(data = as.data.frame(PCA[, 1:3]), 
                                      x = ~PC1, 
                                      y = ~PC2, 
                                      z = ~PC3, 
                                      text = texts, 
                                      color = as.character(colors), 
                                      marker = list(size = 2)))
 
}

spaces <- function(depth) {
  strrep("  ", depth)
}

writeVector <- function(fd, depth, v, title) {
  cat(spaces(depth), title, ": ", file = fd, sep = "")
  cat("[", paste0(v, collapse = ", "), "]\n", file = fd, sep = "")
}

writeMatrix <- function(fd, depth, m, title) {
  cat(spaces(depth), title, ":\n", sep = "", file = fd)
  ss <- spaces(depth + 1)
  rowF <- function(x) paste0(ss, "- [", paste0(x, collapse = ", "), "]", "\n")
  cat(paste0(apply(m, 1, rowF)), file = fd, sep = "")
}

writeCluster <- function(fd, cluster, i) {
  cat("  - cluster: ", i, "\n", file = fd)
  writeVector(fd, 2, cluster$mean, 'mean')
  writeMatrix(fd, 2, cluster$US, "US")
  writeMatrix(fd, 2, cluster$counts, "counts")
}

writeYaml <- function(clusterResults, clustering, variants, 
                      outputFileName, title) {
  if(missing(outputFileName)){
    stop("Output file name missing")
  }
  
  fd <- file(outputFileName, open = "w")
  on.exit(close(fd))
  
  write <- function(...) {
    cat(..., file = fd, sep = "")
  }
  
  write("title: ", title, "\n")
  write("version: ", toString(utils::packageVersion(utils::packageName())), "\n")
  write("salt: ", stringi::stri_rand_strings(1, 16, "[0-9a-f]"), "\n")
  titleF <- function(x) x[["title"]]
  writePopulationStructure(clustering$hier, clustering$classes, fd)
  write("\nvariants:\n")
  for (v in variants) {
    write("  - ", v, "\n")
  }
  write("population:\n")
  for (i in 1:length(clusterResults)) {
    writeCluster(fd, clusterResults[[i]], i)
  }
}


drop <- function(pca, knn_rate, mvn_rate) {
  n <- ncol(pca)
  knn_n <- ceiling(n * knn_rate)
  mvn_n <- ceiling(n * mvn_rate)
  if (knn_n + mvn_n > n - nrow(pca)) {
    stop("Drop rates are too high.")
  }
  
  ids <- 1:n
  
  knn_drop <- adamethods::do_knno(t(pca), 5, knn_n)  
  pca <- pca[, -knn_drop]
  ids <- ids[-knn_drop]
  
  ids[normal_subsample(pca, n - knn_n - mvn_n)$points]
  
}

#' Prepare instance and write yaml file from QC-ed gmatrix.
#' 
#' @param gmatrix gmatrix with imputed missing values.
#' @param imputationResults matrix of the same dimensions as \code{gmatrix} 
#' filled with values TRUE/FALSE indicating which genotypes were imputed.
#' @param controlsU Reference U matrix from SVD decomposition of centered
#' genotype matrix of control dataset or joint dataset.
#' @param meanControl numeric vector representing mean per-variant genotype, 
#' that were taken as a center while costructing \code{controlsU} matrix.
#' @param outputFileName name of the YAML file to output.
#' @param clusters clustering object.
#' @param title title to put as a first line in yaml file.
#' @param MAC minor allele count used for QC.
#' @param MAF minor allele frequency used for QC.
#' @param knn_drop fraction of samples that would be excluded from dataset at
#' each level of hierarchy clustering by applying knn anomaly detection algorithm.
#' @param normalize_drop fraction of samples to be excluded from dataset at
#' each level of hierarchy clustering by applying function \code{normal_subsample} 
#' from the package.
#' @import mclust
#' @export
prepareInstance <- function (gmatrix, imputationResults, controlsU, meanControl, 
                             outputFileName, clusters = NULL, 
                             title = "DNAScoreInput", MAC = 10, MAF = 0.01, 
                             knn_drop = 0.05, normalize_drop = 0.05) {
  if (is.null(clusters)) {
    classes <- rep("Main", ncol(gmatrix))
    clusters <- clustering(setNames(classes, colnames(gmatrix)), 
                           "Main")
  }
  stopifnot("clustering" %in% class(clusters))
  gmatrixForCounts <- gmatrix
  gmatrixForCounts[which(imputationResults, arr.ind = TRUE)] <- NA
  caseCounts <- list()
  clusters$hier$Do(function(node) {
    i <- ifelse(is.null(node$id), return(), node$id)
    cluster_leaves <- sapply(node$leaves, function(x) x$id)
    cluster <- which(clusters$samples %in% cluster_leaves)
    caseCounts[[i]] <<- genotypesToCounts(gmatrixForCounts[, cluster])
  })
  passVariants <- lapply(caseCounts, checkAlleleCounts, mac = MAC, maf = MAF)
  passVariants <- lapply(passVariants, which)
  passVariants <- Reduce(intersect, passVariants)
  if (length(passVariants) < 100) {
    stop("Less than 100 variants left after QC.")
  }
  if (length(passVariants) < 500) {
    warning("Less than 500 variants left after QC.")
  }
  gmatrix <- gmatrix[passVariants, ]
  gmatrixForCounts <- gmatrixForCounts[passVariants, ]
  
  names(meanControl) <- rownames(controlsU)
  controlsU <- controlsU[rownames(gmatrix), ]
  meanControl <- meanControl[rownames(gmatrix)]
  
  numberOfClusters <- 2 * length(clusters$classes) - 1
  clusterResults <- vector("list", numberOfClusters)
  
  cluster_ids <- list()
  gm_names <- c()
  if (is.null(colnames(gmatrix))) {
    gm_names <- 1:ncol(gmatrix)
  } else {
    gm_names <- colnames(gmatrix)
  }
  
  clusters$hier$Do(function(node) {
    i <- ifelse(is.null(node$id), return(), node$id)
    cluster_leaves <- sapply(node$leaves, function(x) x$id)
    
    initial_cluster <- which(clusters$samples %in% cluster_leaves)
    clusterGenotypes <- t(controlsU) %*% (gmatrix[, initial_cluster] - meanControl)
    
    cluster <- initial_cluster[drop(clusterGenotypes, knn_drop, normalize_drop)]
    clusterGenotypes <- t(controlsU) %*% (gmatrix[, cluster] - meanControl)
    
    clusterGenotypesForCounts <- gmatrixForCounts[, cluster]
    clusterMeans <- rowMeans(clusterGenotypes)
    k <- min(nrow(clusterGenotypes), ncol(clusterGenotypes))
    clusterGenotypes <- clusterGenotypes - clusterMeans
    svdResult <- suppressWarnings(RSpectra::svds(A = clusterGenotypes, k = k))
    US <- svdResult$u %*% diag(svdResult$d)
    US <- US / sqrt(length(cluster))
    
    counts <- genotypesToCounts(clusterGenotypesForCounts)
    cluster_ids[[i]] <<- gm_names[cluster]
    clusterResults[[i]] <<- list(US = US, counts = counts, 
                                mean = clusterMeans, 
                                title = clusters$classes[i])
    save(i, cluster_leaves, initial_cluster, cluster, file = paste0(node$id, ".txt"))
  })
  writeYaml(clusterResults, clusters, variants = rownames(gmatrix), 
            outputFileName = outputFileName, title = title)
  cluster_ids
}

processHierarchy <- function(hier) {
  if (names(hier) == "split") {
    list(left = processHierarchy(hier$split$left), 
         right = processHierarchy(hier$split$right), 
         type = "split", id = hier$split$id)
  } else {
    c(hier$cluster, list(type = "leaf"))
  }
}

getNamesFromHier <- function(hier) {
  if (!(names(hier) %in% c("split", "cluster"))) {
    userError("Incorrect hierarchy section in the input YML file.")
  }
  if (names(hier) == "split") {
    if (length(hier$split) != 3) {
      userError("Incorrect hierarchy section in the input YML file.")
    }
    ret <- c(getNamesFromHier(hier$split$left), getNamesFromHier(hier$split$right))
    name <- paste(ret, collapse = " & ")
    c(ret, setNames(name, hier$split$id))
  }  else {
    if (!setequal(names(hier$cluster), c("id", "name"))) {
      userError("Incorrect hierarchy section in the input YML file.")
    }
    stats::setNames(hier$cluster$name, hier$cluster$id)
  } 
}

#' Read instance from YML file.
#' 
#' @param filename YML file.
#' @return a list of four: title of dataset, version, list of variants, hierarchy and
#' population description.
#' @export
readInstanceFromYml <- function(filename) {
  tryCatch(inst <- yaml::read_yaml(filename), 
           error = function(e) userError("File is not a correct YML file"))
  if (!setequal(names(inst), c("title", "salt", "version", "hierarchy", "variants", "population"))) {
    userError("Incorrect input YML file: missing or extra sections.")
  }
  i <- 1
  inst$population <- lapply(inst$population, function(x) {
    err <- function(x) {
      userError(paste0(x, " The error occurred while processing population #", i, "."))
    }
    for (obj in c("counts", "US")) {
      raw <- x[[obj]]
      if (length(raw) == 0) {
        err(paste("Missing", obj, "matrix."))
      }
      x[[obj]] <- do.call(rbind, x[[obj]])
      mode(x[[obj]]) <- "numeric"
      if (any(is.na(x[[obj]]))) {
        err("US and counts matrices in input must contain only numeric values.")
      }
      if (obj == "counts" && nrow(x[[obj]]) != length(inst$variants)) {
        err("Number of rows in counts matrix must match number of variants 
                   in the input YML file.")
      }
      if (length(unique(sapply(raw, length))) != 1) {
        err(paste0("Malformed one of the ", obj, " matrices in the input YML file."))
      }
    }
    if (ncol(x$counts) != 3) {
      userError("Number of columns of counts matrix in the YML file must be three.")
    }
    x$id <- i
    i <<- i + 1
    x
  })
  if (any(sapply(inst$population, function(x) x$cluster) != 1:length(inst$population))) {
    userError("List of populations must go in natural order by their ids")
  }
  inst$names <- getNamesFromHier(inst$hierarchy)
  if (!setequal(names(inst$names), 1:length(inst$population))) {
    userError("Hierarchy section in input YML file have entries with unexpected IDs")
  }
  inst$hierarchy <- processHierarchy(inst$hierarchy)
  inst
}