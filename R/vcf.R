transformToTabixRegions <- function(x, pattern, outputPattern){
  captures <- regexpr(pattern, x, perl = TRUE)
  if (captures == -1) {
    stop(paste0("Region ", x, " is not correctly formatted"))
  }
  starts <- attr(captures, "capture.start")
  stops <- starts + attr(captures, "capture.length") - 1
  substrings <- sapply(mapply(c, starts, stops, SIMPLIFY = FALSE), 
    function(y) {
      substr(x, y[1], y[2])
    })
  paste0(substrings[outputPattern[1]], ":", substrings[outputPattern[2]], "-",
         substrings[outputPattern[3]])
}

createVCFFromTabixIndex <- function(vcf, variants, regions) {
  variantsPositions <- NULL
  callRatePosition <- NULL
  header <- seqminer::tabix.read.header(vcf)
  t <- tempfile("regions", fileext = ".vcf.gz")
  con <- gzfile(t, "w")
  write(header$header, con)
  if (!is.null(regions)) {
    regionsPattern <- "^[\\s]*chr(\\d{1,2}|X|Y) [\\s]*([\\d]+)[\\s]*([\\d]+)[\\s]*$"
    crt <- function(x) transformToTabixRegions(x, regionsPattern, c(1, 2, 3))
    callRatePositions <- sapply(regions, crt)
    for (region in callRatePositions) {
      write(seqminer::tabix.read(vcf, region), file = con, append = TRUE)
    }
  }
  if (!is.null(variants)) {
    variantsPattern <- "^[\\s]*chr(\\d{1,2}|X|Y)[\\s]*:[\\s]*([\\d]+)[\\s]"
    vart <- function(x) transformToTabixRegions(x, variantsPattern, c(1, 2, 2))
    variantsPositions <- sapply(variants, vart)
    for (var in variantsPositions) {
      write(seqminer::tabix.read(vcf, var), file = con, append = TRUE)
    }
  }
  close(con)
  t
}

#' Scan VCF files 
#' 
#' Scan .vcf or .vcf.gz files in matrix and return genotype matrix, call rate
#' with applied filters. In addition, the method can generate binary and 
#' metadata files for faster access to genotype data.
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
#' (i.e. chr23:1532 T GT). In case of deletion ALT must be "*". If present
#' corrspoding genotype matrix will be returned
#' @param returnGenotypeMatrix logical: whether or not return genotype matrix
#' @param regions the set of regions in format "chr# startPos endPos". For each
#' region call rate will be calculated and corresponding matrix will be returned. 
#' @param binaryPathPrefix the path prefix for binary file prefix_bin and 
#' metadata file prefix_meta. If not NULL corresponding files will be generated.
#' @return list containing genotype matrix and/or call rate matrix if 
#' requested.
#' @export
scanVCF <- function(vcf, DP = 10L, GQ = 20L, samples = NULL,
                    bannedPositions = NULL, variants = NULL, 
                    returnGenotypeMatrix = TRUE, 
                    regions = NULL, binaryPathPrefix = NULL) {
  stopifnot(length(DP) > 0)
  stopifnot(length(GQ) > 0)
  DP <- as.integer(DP)
  GQ <- as.integer(GQ)
  stopifnot(!is.na(DP[0]))
  stopifnot(!is.na(GQ[0]))
  
  stopifnot(file.exists(vcf))
  
  tbi <- paste0(vcf, ".tbi")
  if (is.null(binaryPathPrefix) && file.exists(tbi) && 
      !(length(variants) == 0) && returnGenotypeMatrix) {
    vcf <- createVCFFromTabixIndex(vcf, variants, regions)
  } else {
    tbi <- NULL
  }
  
  fixChar <- function(x) if(is.null(x)) character(0) else x
  samples <- fixChar(samples)
  bannedPositions <- fixChar(bannedPositions)
  variants <- fixChar(variants)
  regions <- fixChar(regions)
  binaryPathPrefix <- fixChar(binaryPathPrefix)
  
  tryCatch(
    res <- parse_vcf(vcf, samples, bannedPositions, variants, DP, GQ, 
                   returnGenotypeMatrix, regions, binaryPathPrefix),
    error = function(c) {
      suffix <- ""
      if (!is.null(tbi)) {
        suffix <- paste0("Note: since .tbi is provided a temp file containing ",
                "only necessary regions has been created(", vcf, "). See ",
                "the temp file for more information about the error.\n")
      }
      stop(paste0(conditionMessage(c), "\n", suffix))
  })
  
  if (!is.null(res$genotype)) {
      colnames(res$genotype) <- res$samples
  }
  if (!is.null(res$callrate)) {
      colnames(res$callrate) <- res$samples
      rownames(res$callrate) <- regions
  }
  res
}

#' Scan binary variant file
#' 
#' Scan binary file produced by \code{\link{scanVCF}} and collect required 
#' per variant allele counts statistics after applying quality control filters.
#' @inheritParams scanVCF
#' @param binaryFile the name of binary file
#' @param metafile the name of metadata file
#' @param samples the set of samples to be analyzed
#' @param MAC minor allele count
#' @param MAF minor allele frequency
#' @return data.frame with four columns: variant(position, reference allele, 
#' alternative allele); number of samples with 
#' both reference, one reference and alternative allele indicated in first column,
#' and both alternative respectively.
#' @export
scanBinaryFile <- function(binaryFile, metafile, samples, variants, DP = 10, GQ = 20,
                           MAC = 1, MAF = 0.00) {
  stopifnot(file.exists(binaryFile))
  stopifnot(file.exists(metafile))
  res <- parse_binary_file(variants, samples, binaryFile, metafile, DP, GQ);
  total <- as.integer(res["total"])
  res["total"] <- NULL
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  actualMAC <- apply(res[, 2:4], 1, function(x) 2 * min(x[1], x[3]) + x[2])
  res <- res[actualMAC >= MAC & actualMAC / 2 * nrow(res) >= MAF, ]
  res <- res[res$HOM_REF + res$HET + res$HOM >= 0.9 * total, ]
  rownames(res) <- NULL
  res
}
