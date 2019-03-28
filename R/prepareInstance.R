#' Prepare an instance to be submitted to website dnascore.net
#' 
#' The method generates all required files
#' @param vcf vcf file name
#' @param dataset one of predefined in this package datasets. List of supported 
#' datasets: \code{\link{finSwedDataset}}, \code{\link{publicExomesDataset}}.
#' @param outputPathPrefix prefix for ouput files
#' @param vectors number of vectors to be generated
#' @param ... arguments that will be passed to \code{\link{scanVCF}}
#' @param verbose logical
#' @export
prepareInstance <- function(vcf, dataset, outputPathPrefix = "dnascore", 
                            vectors = 10, verbose = FALSE, ...) {
  if (verbose) {
    cat("Parsing genotype matrix...\n")
  }
  genotype <- genotypeMatrixVCF(vcf, variants = dataset$ancestryVariants, 
                                verbose = verbose, ...)
  genotype <- genotype[unique(rownames(genotype)),]
  gmatrix <- replaceMissing(as.matrix(genotype))
  if (nrow(gmatrix) < 2 | ncol(gmatrix) < 2) {
    stop(paste0("Too small genotype matrix. Dimensons: ", nrow(gmatrix), 
                "x", ncol(gmatrix), ". Check DP and GQ filters and ensure ",
                "that vcf file covers required ancestry variants passed as", 
                " dataset argument. Also check quality of vcf file. See", 
                "?scanVCF for more information about default filters."))
  }
  vectorsGM <- min(vectors, ncol(gmatrix) - 1, nrow(gmatrix) - 1)
  if (verbose) {
    cat("Running SVD on genotype matrix...\n")
  }
  svdOut <- RSpectra::svds(gmatrix - rowMeans(gmatrix), k = vectorsGM)
  u <- svdOut$u
  gmatrix <- NULL
  utils::write.table(u, file = paste0(outputPathPrefix, "U.tsv"), 
                     row.names = FALSE, col.names = FALSE)
  u <- NULL
  caseCounts <- genotypesToCounts(genotype)
  colnames(caseCounts)[1] <- "POS\tREF\tALT"
  data.table::fwrite(caseCounts, file = paste0(outputPathPrefix, "CaseCounts.tsv"),
                     quote = FALSE, sep = "\t")
  genotype <- NULL
  caseCounts <- NULL
  
  if (verbose) {
    cat("Parsing VCF file for call rate matrix...\n")
  }
  callrate <- callRateMatrixVCF(vcf, dataset$intervals, ...)
  if (nrow(callrate) == 0) {
    warning("Callrate matrix is empty. Check DP and GQ filters and ensure", 
            "that vcf file contains variants that cover required intervals ", 
            "(see dataset object for intervals). Also check quality of the VCF",
            "file. See ?scanVCF for more information about default filters.")
  } else {
    callrateU <- svd(as.matrix(callrate), nu = min(vectors, nrow(callrate)))$u
  }
}
