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


#' Estimate clusters within case data
#' Clusters are estimated from PC loadings using
#' Gaussian mixed model fitting from Mclust package
#' using ellipsoidal, varying volume, shape, and orientation
#' model ("VVV")
#' @param PCA results of gmatrix principal component analysis
#' @param plotBIC logical whether a plot of Bayesian Information 
#' Criterion vs cluster number should be returned
#' @param plotDendrogram whether dendrogram should be plotted
#' @export
estimateCaseClusters <- function(PCA, plotBIC = FALSE, plotDendrogram = FALSE){
  stopifnot(class(PCA) %in% c("matrix", "data.frame"))
  
  clResults <- mclust::Mclust(data = PCA, G = 1:6, modelNames = "VVV")
  
  if (plotBIC){
    mclust::plot(x = clResults, 
                 what = "BIC", 
                 xlab = "Number of Components", 
                 ylab = "Bayesian Information Criterion")
  }
  
  clCombi <- mclust::clustCombi(clResults)
  collapsing <- dendrogramEstimate(clCombi$combiM)
  
  if (plotDendrogram){
    mclust::combiTree(clCombi, type = "rectangle", yaxis = "entropy")
  }
  
  clResults$classification
}


#' Estimate PCA from gmatrix without missing values 
#' 
#' PCA is estimated using SVD of the mean centered genotype matrix
#' PCA coordinates will be needed to perform clustering of the 
#' case samples to ensure relatively homogenious matching.
#' @param gmatrix Genotype matrix wihout missing values
#' @export
gmatrixPCA <- function(gmatrix, components = 30){
  pca <- RSpectra::svds(A = gmatrix - rowMeans(gmatrix), 
                        k = min(components, ncol(gmatrix)))$v
  
  colnames(pca) <- c(paste("PC", c(1:min(components, ncol(gmatrix))), sep = ""))
  rownames(pca) <- colnames(gmatrix)
  pca
}

#' Filter gmatrix by minimal sample call rate and minimal variant call rate
#' @export
filterGmatrix <- function(gmatrix, imputationResults) {
  nSamples <- ncol(gmatrix)
  nVariants <- nrow(gmatrix)
  sampleCallRate <- apply(imputationResults, MARGIN = 2, function(x){
    length(x[! x]) / length(x)
  })
  gmatrix <- gmatrix[, which(sampleCallRate > minSampleCallRate)]
  
  variantCallRate <- apply(imputationResults, MARGIN = 1, function(x) {
    length(x[!x]) / length(x) 
  })
  
  gmatrix <- gmatrix[which(variantCallRate > minVariantCallRate), ]
  imputationResults <- imputationResults[variantCallRate > minVariantCallRate, 
                                         sampleCallRate > minSampleCallRate]
  
  message("Kept ", ncol(gmatrix), " out of ", nSamples, " individuals")
  message("Kept ", nrow(gmatrix), " out of ", nVariants, " variants")
  list(gmatrix = gmatrix, imputationResults = imputationResults)
}

#' Plot 3d PCA plot 
#' @param PCA
#' @param clusters character vector of cluster assignment for every
#' sample in gmatrix
#' @inheritParams estimateCaseClusters
#' @export 
plotPCA <- function(PCA, clusters = NULL) {
  colors <- if (is.null(clusters)) 1 else clusters 
  texts <- paste("Sample:", rownames(PCA))
  if (!is.null(clusters)) {
    texts <- paste0(texts, "\nCluster: ", clusters)
  }
  plotly::add_markers(plotly::plot_ly(x = PCA[, "PC1"], 
                                      y = PCA[, "PC2"],
                                      z = PCA[, "PC3"], 
                                      text = texts, 
                                      color = colors, 
                                      marker = list(size = 2)))
 
}

spaces <- function(depth) {
  strrep("  ", depth)
}

writeMatrix <- function(fd, depth, m, title) {
  cat(spaces(depth), title, ":\n", sep = "", file = fd)
  ss <- spaces(depth + 1)
  rowF <- function(x) paste0(ss, "- [", paste0(x, collapse = ", "), "]", "\n")
  cat(paste0(apply(m, 1, rowF)), file = fd, sep = "")
}

writeCluster <- function(fd, cluster) {
  cat("  - cluster: ", cluster$title, "\n", file = fd)
  writeMatrix(fd, 2, cluster$US, "US")
  writeMatrix(fd, 2, cluster$counts, "counts")
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

writeYaml <- function(clustering, variants, outputFileName, title) {
  if (length(clusterResults) > 1 && is.null(collapsing) ){
    stop("More than 1 cluster detected, but no collapsing scheme supplied")
  }
  
  if(missing(outputFileName)){
    stop("Output file name missing")
  }
  
  tree <- collapsingToTree(collapsing)
  fd <- file(outputFileName, open = "w")
  on.exit(close(fd))
  cat("title: ", title, "\n", file = fd, sep = "")
  cat("hierarchy:\n  ")
  titleF <- function(x) x[["title"]]
  writePopulationStructure(tree, fd, 1, sapply(clusterResults, titleF))
  cat("\nvariants:\n", file = fd)
  for (v in variants) {
    cat("  - ", v, "\n", sep = "", file = fd)
  }
  cat("population:\n", file = fd, sep = "")
  for (i in 1:length(clusterResults)) {
    writeCluster(fd, clusterResults[[i]])
  }
}

#' prepare instance and write yaml file from QC-ed gmatrix
#' @param gmatrix gmatrix with imputed missing values
#' @param imputationResults matrix indicating which genotypes
#' were imputed
#' @param minVariantCallRate minimal variant call rate
#' @param maxVectors maximum number of principal directions for each cluster
#' to be computed and written to a file
#' @param minSampleCallRate minimal sample call rate
#' @param outputFileName name of the YAML file to output
#' @param title title to put as a first line in yaml file
#' @import mclust
#' @export
prepareInstance <- function(gmatrix, imputationResults, 
                            outputFileName, 
                            clusters = NULL, 
                            maxVectors = 50, 
                            title = "DNAScoreInput"){
  if (is.null(clusters)) {
    classes <- rep("Main", nrow(gmatrix))
    clusters <- clustering(setNames(classes, rownames(gmatrix)), "Main")
  }
  stopifnot("clsutering" %in% class(clusters))
  numberOfClusters <- length(clusters$classes)

  gmatrixForCounts <- gmatrix
  gmatrixForCounts[which(imputationResults, arr.ind = TRUE)] <- NA
  
  caseCounts <- list()
  for (i in 1:numberOfClusters){
    cluster <- which(clusters$classification == i)
    caseCounts[[i]] <- genotypesToCounts(gmatrixForCounts[, cluster])
  }
  
  passVariants <- lapply(caseCounts, checkAlleleCounts)
  passVariants <- lapply(passVariants, which)
  passVariants <- Reduce(intersect, passVariants)
    
  if(length(passVariants) < 2){
    stop("Less than 2 good variants were selected")
  }
  
  gmatrix <- gmatrix[passVariants, ]
  clusterResults <- vector("list", numberOfClusters)
  for (i in 1:numberOfClusters){
    cluster <- which(clusters$classification == i) 
    clusterGenotypes <- gmatrix[, cluster]
    k <- min(length(cluster), nrow(gmatrix), maxVectors)
    svdResult <- suppressWarnings(RSpectra::svds(A = clusterGenotypes, 
                                                 k = k))
    US <- svdResult$u %*% diag(svdResult$d)
    counts <- genotypesToCounts(clusterGenotypes)
    clusterResults[[i]] <- list(US = US[, -1], counts = counts, 
                                title = paste0("cluster", i))
  }
  
  writeYaml(clusterResults = clusterResults, 
            variants = rownames(gmatrix),
            outputFileName = outputFileName, 
            clusters = clusters, title = title)
  
}
