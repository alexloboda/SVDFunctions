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
  
  clResults <- mclust::Mclust(data = PCA, G = 1:20, modelNames = "VVV")
  
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
  
  list(clResults = clResults, collapsing = collapsing)
}


#' Estimate PCA from gmatrix without missing values and make 3dplot
#' 
#' PCA is estimated using SVD of the mean centered genotype matrix
#' PCA coordinates will be needed to perform clustering of the 
#' case samples to ensure relatively homogenious matching.
#' @param gmatrix Genotype matrix wihout missing values
#' @param clusters character vector of cluster assignment for every
#' sample in gmatrix
#' @import plotly
#' @export
gmatrixPCA <- function(gmatrix, clusters = NULL){
  pca <- RSpectra::svds(A = gmatrix - rowMeans(gmatrix), 
                        k = min(30, ncol(gmatrix)))$v
  
  colnames(pca) <- c(paste("PC", c(1:min(30, ncol(gmatrix))), sep = ""))
  colors <- if (is.null(clusters)) 1 else clusters 
  texts <- paste("Sample:", colnames(gmatrix))
  p <- plotly::add_markers(plotly::plot_ly(x = pca[, "PC1"], 
                                           y = pca[, "PC2"],
                                           z = pca[, "PC3"], 
                                           text = texts, 
                                           color = colors, 
                                           marker = list(size = 2)))
  
  list(PCA = pca, plot = p)  
}

spaces <- function(depth) {
  strrep("  ", depth)
}

cats <- function(file, depth, ...) {
  cat(spaces(depth), ..., file = file, sep = "")
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

writePopulationStructure <- function(tree, fd, depth, titles) {
  if (length(tree) == 1) {
    cats(fd, 0, "cluster:\n")
    id <- tree[[1]]
    cats(fd, depth + 1, "id: ", id, "\n")
    cats(fd, depth + 1, "name: ", titles[id], "\n")
  } else {
    cats(fd, 0, "split:\n")
    for (i in 1:length(tree)) {
      cats(fd, depth + 1, "- ")
      writePopulationStructure(tree[[i]], fd, depth + 2, titles)
    }
  }
}

writeYaml <- function(clusterResults, variants, outputFileName, collapsing, 
                      title) {
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
prepareInstance <- function(gmatrix, imputationResults, outputFileName, 
                            maxVectors = 50, 
                            minVariantCallRate = 0.95,
                            minSampleCallRate = 0.95,
                            title = "DNAScoreInput"){
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
  
  pca <- RSpectra::svds(A = gmatrix - rowMeans(gmatrix), 
                        k = min(30, ncol(gmatrix), nrow(gmatrix)))$v
  
  clResults <- mclust::Mclust(data = pca, G = 1:20, modelNames = "VVV")
  numberOfClusters <- length(unique(clResults$classification))

  gmatrixForCounts <- gmatrix
  gmatrixForCounts[which(imputationResults, arr.ind = TRUE)] <- NA
  if (numberOfClusters > 1){
    clCombi <- mclust::clustCombi(clResults)
    collapsing <- dendrogramEstimate(clCombi$combiM)
    
    case_counts <- vector("list", numberOfClusters)
    for (i in 1:numberOfClusters){
      cluster <- which(clResults$classification == i)
      case_counts[[i]] <- genotypesToCounts(gmatrixForCounts[, cluster])
    }
    
    pass_variants <- lapply(case_counts, checkAlleleCounts)
    pass_variants <- lapply(pass_variants, which)
    pass_variants <- Reduce(intersect, pass_variants)
  } else{
    case_counts <- genotypesToCounts(gmatrixForCounts) 
    pass_variants <- which(checkAlleleCounts(case_counts))
  }
  
  if(length(pass_variants) < 2){
    stop("Less than 2 good variants were selected")
  }
  gmatrix <- gmatrix[pass_variants, ]
  clusterResults <- vector("list", numberOfClusters)
  for (i in 1:numberOfClusters){
    cluster <- which(clResults$classification == i) 
    k <- min(length(cluster), nrow(gmatrix), maxVectors)
    svdResult <- suppressWarnings(RSpectra::svds(A = gmatrix[, cluster], 
                                                 k = k))
    US <- svdResult$u %*% diag(svdResult$d)
    counts <- genotypesToCounts(gmatrix[, cluster])
    clusterResults[[i]] <- list(US = US[, -1], counts = counts, 
                                title = paste0("cluster", i))
  }
  names(clusterResults) <- c(1:numberOfClusters)
  
  if(numberOfClusters == 1){
    collapsing <- NULL
  }
  writeYaml(clusterResults = clusterResults, 
            variants = rownames(gmatrix),
            outputFileName = outputFileName, 
            collapsing = collapsing, title = title)
  
  message("Kept ", ncol(gmatrix), " out of ", nSamples, " individuals")
  message("Kept ", nrow(gmatrix), " out of ", nVariants, " variants")
}
