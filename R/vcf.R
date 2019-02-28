#' Scan VCF files 
#' 
#' Scan .vcf or .vcf.gz files in matrix and return genotype matrix, call rate
#' with applied filters. Also can generate binary and metadata files for
#' faster access to genotype data.
#' @param vcf the name of file to read, can be plain text VCF file as well
#' as compressed with gzip or zlib headers.
#' @param DP integer: minimum required read depth for position to be considered,
#' otherwise assumed as missing.
#' @param GQ integer: minimum required genotype quality for position to be 
#' considered, otherwise assumed as missing.
#' @param samples the set of samples to be scanned and returned
#' @param bannedPositions the set of positions in format "chr#:#" that 
#' must be eliminated from consideration. 
#' @param variants the set of variants in format "chr#:# REF ALT"
#' (i.e. chr23:1532 T GT). In case of deletion ALT must be "*". 
#' @param returnGenotypeMatrix logical: if TRUE genotype matrix will be returned
#' @param regions the set of regions in format "chr# startPos endPos". For each
#' region call rate will be calculated and corresponding matrix will be returned. 
#' @param binaryPathPrefix the path prefix for binary file prefix_bin and 
#' metadata file prefix_meta. If not NULL corresponding files will be generated.
#' @return list containing genotype matrix and/or call rate matrix if 
#' requested.
#' @export
scanVCF <- function(vcf, DP = 10L, GQ = 20L, samples = NULL, 
                    bannedPositions = NULL, variants = NULL, 
                    returnGenotypeMatrix = TRUE, regions = NULL,
                    binaryPathPrefix = NULL) {
  stopifnot(length(DP) > 0)
  stopifnot(length(GQ) > 0)
  DP <- as.integer(DP)
  GQ <- as.integer(GQ)
  stopifnot(!is.na(DP[0]))
  stopifnot(!is.na(GQ[0]))
  
  stopifnot(file.exists(vcf))
  
  fixChar <- function(x) if(is.null(x)) character(0) else x
  samples <- fixChar(samples)
  bannedPositions <- fixChar(bannedPositions)
  variants <- fixChar(variants)
  regions <- fixChar(regions)
  binaryPathPrefix <- fixChar(binaryPathPrefix)
  
  res <- parse_vcf(vcf, samples, bannedPositions, variants, DP, GQ, 
                   regions, returnGenotypeMatrix, binaryPathPrefix);
  
  if (!is.null(res$genotype)) {
      colnames(res$genotype) <- res$samples
  }
  if (!is.null(res$callrate)) {
      colnames(res$callrate) <- res$samples
      rownames(res$callrate) <- regions
  }
  res
}

#' @export
scanBinaryFile <- function(binaryFile, metafile, samples, variants, DP = 10, GQ = 20) {
  stopifnot(file.exists(binaryFile))
  stopifnot(file.exists(metafile))
  res <- parse_binary_file(variants, samples, binaryFile, metafile, DP, GQ);
  as.data.frame(res)
}
