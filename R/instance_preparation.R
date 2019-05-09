#' Estimate dendrogram for case clustering
#' uses combiM output of clustCombi() function
#' from Mclust package
#' @param combiM combiM output of clustCombi() function
#'
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
#' @param plot logical whether a plot of Bayesian Information 
#' Criterion vs cluster number should be returned
#' @export
estimateCaseClusters <- function(PCA, plotBIC = FALSE, plotDendrogram = FALSE){
  stopifnot(class(PCA) %in% c("matrix", "data.frame"))
  
  clResults <- mclust::Mclust(data = PCA, G = 1:20, modelNames = "VVV")
  
  if (plotBIC){
    mclust::plot.Mclust(x = clResults, 
                        what = "BIC", 
                        xlab = "Number of Components", 
                        ylab = "Bayesian Information Criterion")
  }
  
  clCombi <- mclust::clustCombi(clResults)
  collapsing <- dendrogramEstimate(clCombi)
  
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
#' @export
gmatrixPCA <- function(gmatrix, clusters = NULL){
  pca <- RSpectra::svds(A = gmatrix - rowMeans(gmatrix), k = min(30, ncol(gmatrix)))$v
  pca <- cbind(pca, clusters)
  colnames(pca) <- c(paste("PC", c(1:min(30, ncol(gmatrix))), sep = ""),
                     "Cluster")
  p <- plotly::add_markers(plotly::plot_ly(x = PC1, 
                                           y = PC2,
                                           z = PC3, 
                                           color = Cluster,
                                           marker = list(size = 2 )))
  
  list(PCA = pca, plot = p)  
}


#' QC gmatrix
#' performs callrate estimate on samples and variants
#' checks for Hardy-Weinberg equilibrium and finds consistent
#' set of variants across all clusters
#' 

gmatrixQC <- function(gmatrix, imputationResults, clusters, minVariantCallRate = 0.95,
                      minSampleCallRate = 0.95){
  sampleCallRate <- apply(imputationResults, MARGIN = 2, function(x){
    length(x[x == 0]) / length(x) })
  gmatrix <- gmatrix[, which(sampleCallRate > minSampleCallRate)]
  
  variantCallRate <- apply(imputationResults, MARGIN = 1, function(x){
    length(x[x == 0]) / length(x) })
  gmatrix <- gmatrix[which(variantCallRate > minVariantCallRate), ]
  
  pass_variants <- vector("list", length(unique(clusters)))
  for (i in 1:length(unique(clusters))){
    pass_variants <- genotypesToCounts()
  }
  
}

#' Write yaml file given the list of US and case counts for 
#' all clusters
#' @param clusterResults list of 2 objects - list of US and list of 
#' case counts list(US, case_counts)
#' @param outputFileName name of the yaml file to write to
#' @param collapsing matrix representing cluster collapsing
#' scheme
#' @param title dataset name

writeYaml <- function(clusterResults, variants, outputFileName, collapsing, title = "DNAscoreInput"){
  if (length(clusterResults) > 1 & missing(collapsing)){
    stop("More than 1 cluster detected, but no collapsing scheme supplied")
  }
  if(missing(outputFileName)){
    stop("Output file name missing")
  }
  for (i in 1:length(clusterResults)){
    clusterResults[[i]] <- data.frame(US = paste("[", 
                                                 apply(clusterResults[[i]]$US,
                                                       MARGIN = 1,
                                                       paste, collapse = ", "),
                                                 "]", sep = ""),
                                      Counts = paste("[", 
                                                 apply(clusterResults[[i]]$counts,
                                                       MARGIN = 1,
                                                       paste, collapse = ", "),
                                                 "]", sep = ""))
    # clusterResults[[i]][[1]] <- 
    # data.frame(US = paste("[", 
    #                       apply(clusterResults[[i]]$US, 
    #                             MARGIN = 1, 
    #                             paste, collapse = ", "), 
    #                       "]", sep = ""))
    # clusterResults[[i]][[2]] <- 
    #   data.frame(Counts = paste("[", 
    #                       apply(clusterResults[[i]]$counts, 
    #                             MARGIN = 1, 
    #                             paste, collapse = ", "), 
    #                       "]", sep = ""))
    # names(clusterResults[[i]]) <- ""
  }
  if(missing(collapsing)){
    collapsingScheme <- clusterResults
  } else{
    collapsingScheme <- vector("list", nrow(collapsing))
    for (i in 1:length(collapsingScheme)){
      collapsingScheme[[i]] <- vector("list", 2)
      if (collapsing[i, 1] < 0){
        collapsingScheme[[i]][1] <- clustersResults[[-collapsing[i, 1]]]
      } else{
        collapsingScheme[[i]][1] <- collapsingScheme[[collapsing[i, 1]]]
      }
      if (collapsing[i, 2] < 0){
        collapsingScheme[[i]][2] <- clustersResults[[-collapsing[i, 2]]]
      } else{
        collapsingScheme[[i]][2] <- collapsingScheme[[collapsing[i, 2]]]
      }
    }
  }
  
  collapsingScheme
  yaml::write_yaml(list(Population = collapsingScheme,
                        Variants = variants), 
                   file = outputFileName,
                   column.major = TRUE)
}


#' prepare instance and write yaml file from QC-ed gmatrix
#' @param gmatrix gmatrix with imputed missing values
#' @param imputationResults matrix indicating which genotypes
#' were imputed
#' @param minVariantCallRate minimal variant call rate
#' @param minSampleCallRate minimal sample call rate
#' @import mclust
#' @export

prepareInstance1 <- function(gmatrix, imputationResults, minVariantCallRate = 0.95,
                             minSampleCallRate = 0.95, outputFileName,
                             title = "DNAScoreInput"){
  nSamples <- ncol(gmatrix)
  nVariants <- nrow(gmatrix)
  sampleCallRate <- apply(imputationResults, MARGIN = 2, function(x){
    length(x[! x]) / length(x) })
  gmatrix <- gmatrix[, which(sampleCallRate > minSampleCallRate)]
  
  
  variantCallRate <- apply(imputationResults, MARGIN = 1, function(x){
    length(x[! x]) / length(x) })
  gmatrix <- gmatrix[which(variantCallRate > minVariantCallRate), ]
  imputationResults <- 
    imputationResults[which(variantCallRate > minVariantCallRate), 
                      which(sampleCallRate > minSampleCallRate)]
  
  pca <- RSpectra::svds(A = gmatrix - rowMeans(gmatrix), 
                        k = min(30, ncol(gmatrix), nrow(gmatrix)))$v
  
  clResults <- mclust::Mclust(data = pca, G = 1:20, modelNames = "VVV")
  numberOfClusters <- length(unique(clResults$classification))

  gmatrixForCounts <- gmatrix
  gmatrixForCounts[which(imputationResults, arr.ind = TRUE)] <- NA
  if (numberOfClusters > 1){
    clCombi <- mclust::clustCombi(clResults)
    collapsing <- dendrogramEstimate(clCombi)
    
    case_counts <- vector("list", numberOfClusters)
    for (i in 1:numberOfClusters){
      case_counts[[i]] <- 
        genotypesToCounts(gmatrixForCounts[, which(clResults$classification == i)])
    }
    pass_variants <- lapply(case_counts, checkAlleleCounts)
    pass_variants <- lapply(pass_variants, function(x){
      which(x)
    })
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
  if (numberOfClusters > 1) {
    for (i in 1:numberOfClusters){
      
      svdResult <- 
        RSpectra::svds(A = gmatrix[, which(clResults$classification == i)], 
                       k = min(ncol(gmatrix[, which(clResults$classification == i)]),
                               nrow(gmatrix[, which(clResults$classification == i)])))
      US <- svdResult$u %*% diag(svdResult$d)
      counts <- genotypesToCounts(gmatrix[, which(clResults$classification == i)])
      clusterResults[[i]] <- list(US = US[, -1], counts = counts)
    }
    names(clusterResults) <- c(1:numberOfClusters)
  } else {
    svdResult <- 
      RSpectra::svds(A = gmatrix, 
                     k = min(ncol(gmatrix),
                             nrow(gmatrix)))
    US <- svdResult$u %*% diag(svdResult$d)
    counts <- genotypesToCounts(gmatrix)
    clusterResults[[i]] <- list(US = US[, -1], counts = counts)
  }
  
  if(numberOfClusters == 1){
    writeYaml(clusterResults = clusterResults,
              variants = rownames(gmatrix),
              outputFileName = outputFileName, 
              title = title)
  } else{
    writeYaml(clusterResults = clusterResults, 
              variants = rownames(gmatrix),
              outputFileName = outputFileName, 
              collapsing = collapsing, title = title)
  }
  
  message("Kept ", ncol(gmatrix), " out of ", nSamples, " individuals")
  message("Kept ", nrow(gmatrix), " out of ", nVariants, " variants")
  #list(population = clusterResults, variants = rownames(gmatrix))
    
}