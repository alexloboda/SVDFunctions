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

#' Estimate clusters within case data
#' Clusters are estimated from PC loadings using
#' Gaussian mixed model fitting from Mclust package
#' using ellipsoidal, varying volume, shape, and orientation
#' model ("VVV")
#' @param PCA results of gmatrix principal component analysis
#' @param plotBIC logical whether a plot of Bayesian Information 
#' Criterion vs cluster number should be returned
#' @param plotDendrogram whether dendrogram should be plotted
#' @param clusters maximum number of clusters
#' @export
estimateCaseClusters <- function(PCA, plotBIC = FALSE, plotDendrogram = FALSE, 
                                 clusters = 20){
  stopifnot(class(PCA) %in% c("matrix", "data.frame"))
  
  clResults <- mclust::Mclust(data = PCA, G = 1:clusters, modelNames = "VVV")
  
  if (plotBIC){
    mclust::plot.Mclust(x = clResults, 
                        what = "BIC", 
                        xlab = "Number of Components", 
                        ylab = "Bayesian Information Criterion")
  }
  
  clusters <- clResults$G
  if (clusters > 1) {
    clCombi <- mclust::clustCombi(clResults)
    collapsing <- dendrogramEstimate(clCombi$combiM)
  } else {
    collapsing <- NULL
  }
  
  if (plotDendrogram){
    mclust::combiTree(clCombi, type = "rectangle", yaxis = "entropy")
  }
  
  clustering(clResults$classification, collapsingToTree(collapsing))
}


#' Estimate PCA from gmatrix without missing values 
#' 
#' PCA is estimated using SVD of the mean centered genotype matrix
#' PCA coordinates will be needed to perform clustering of the 
#' case samples to ensure relatively homogenious matching.
#' @param gmatrix Genotype matrix wihout missing values
#' @param components Number of principal components to be computed
#' @export
gmatrixPCA <- function(gmatrix, components = 10){
  pca <- RSpectra::svds(A = gmatrix - rowMeans(gmatrix), 
                        k = min(components, ncol(gmatrix)))$v
  
  colnames(pca) <- c(paste("PC", c(1:min(components, ncol(gmatrix))), sep = ""))
  rownames(pca) <- colnames(gmatrix)
  pca
}

#' Filter gmatrix by minimal sample call rate and minimal variant call rate
#' @param minVariantCallRate minimal variant call rate
#' @param minSampleCallRate minimal sample call rate
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

#' Plot 3d PCA plot 
#' @param PCA pca of genotype matrix returned by function gmatrixPCA
#' @inheritParams prepareInstance
#' @inheritParams estimateCaseClusters
#' @export 
plotPCA <- function(PCA, clusters = NULL) {
  colors <- if (is.null(clusters)) 1 else clusters$classes[clusters$samples]
  texts <- paste("Sample:", rownames(PCA))
  if (!is.null(clusters)) {
    texts <- paste0(texts, "\nCluster: ", clusters$samples)
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
  write("hierarchy:\n  ")
  titleF <- function(x) x[["title"]]
  writePopulationStructure(clustering$hier, clustering$classes, fd, 1)
  write("\nvariants:\n")
  for (v in variants) {
    write("  - ", v, "\n")
  }
  write("population:\n")
  for (i in 1:length(clusterResults)) {
    writeCluster(fd, clusterResults[[i]])
  }
}

#' prepare instance and write yaml file from QC-ed gmatrix
#' @param gmatrix gmatrix with imputed missing values
#' @param imputationResults matrix indicating which genotypes
#' were imputed
#' @param maxVectors maximum number of principal directions for each cluster
#' to be computed and written to a file
#' @param outputFileName name of the YAML file to output
#' @param title title to put as a first line in yaml file
#' @param clusters clustering object
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
  
  stopifnot("clustering" %in% class(clusters))
  numberOfClusters <- length(clusters$classes)

  gmatrixForCounts <- gmatrix
  gmatrixForCounts[which(imputationResults, arr.ind = TRUE)] <- NA
  
  caseCounts <- list()
  for (i in 1:numberOfClusters){
    cluster <- which(clusters$samples == i)
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
    cluster <- which(clusters$samples == i) 
    clusterGenotypes <- gmatrix[, cluster]
    k <- min(length(cluster), nrow(gmatrix), maxVectors)
    svdResult <- suppressWarnings(RSpectra::svds(A = clusterGenotypes, 
                                                 k = k))
    US <- svdResult$u %*% diag(svdResult$d)
    counts <- genotypesToCounts(clusterGenotypes)
    clusterResults[[i]] <- list(US = US, counts = counts, 
                                title = i)
  }
  
  writeYaml(clusterResults, clusters, 
            variants = rownames(gmatrix),
            outputFileName = outputFileName, 
            title = title)
  
}