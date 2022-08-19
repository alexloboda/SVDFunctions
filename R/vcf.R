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

createVCFFromTabixIndex <- function(vcf, variants, regions, verbose) {
  variantsPositions <- NULL
  callRatePosition <- NULL
  header <- seqminer::tabix.read.header(vcf)
  t <- tempfile("regions", fileext = ".vcf.gz")
  if (verbose) {
    cat(paste("Creating temp VCF file for specific regions:", t, "\n"))
  }
  con <- gzfile(t, "w")
  write(header$header, con)
  if (verbose) {
    pb <- utils::txtProgressBar(max = length(regions) + length(variants), 
                                style = 3)
  }
  cnt <- 0
  if (!is.null(regions)) {
    regionsPattern <- "^[\\s]*chr(\\d{1,2}|X|Y) [\\s]*([\\d]+)[\\s]*([\\d]+)[\\s]*$"
    crt <- function(x) transformToTabixRegions(x, regionsPattern, c(1, 2, 3))
    callRatePositions <- sapply(regions, crt)
    for (region in callRatePositions) {
      write(seqminer::tabix.read(vcf, region), file = con, append = TRUE)
      if (verbose) {
        cnt <- cnt + 1
        utils::setTxtProgressBar(pb, cnt)
      }
    }
  }
  if (!is.null(variants)) {
    variantsPattern <- "^[\\s]*chr(\\d{1,2}|X|Y)[\\s]*:[\\s]*([\\d]+)[\\s]"
    vart <- function(x) transformToTabixRegions(x, variantsPattern, c(1, 2, 2))
    variantsPositions <- sapply(variants, vart)
    for (var in variantsPositions) {
      write(seqminer::tabix.read(vcf, var), file = con, append = TRUE)
      if (verbose) {
        cnt <- cnt + 1
        utils::setTxtProgressBar(pb, cnt)
      }
    }
  }
  if (verbose) {
    close(pb)
  }
  close(con)
  t
}

#' Scan VCF file for call rate matrix.
#' 
#' The method parses VCF file and calculate call rate for all specified regions.
#' @inheritParams scanVCF
#' @return numeric matrix: ratio of non-missing values to all variants in the
#' region presented in VCF file.
#' @examples
#' vcf <- "CEU.exon.2010_09.genotypes.vcf.gz"
#' filepath <- system.file("extdata", vcf, package = "SVDFunctions")
#' 
#' regions <- c("chr1 1108138 36000000", "chr20 33521213 61665700") 
#' samples <- c("NA06989", "NA10847", "NA11840", "NA12873")
#' 
#' callRateMatrixVCF(filepath, regions, GQ = 0, samples = samples)
#' @export
callRateMatrixVCF <- function(vcf, regions, DP = 10L, GQ = 20L, samples = NULL, 
                  bannedPositions = NULL, verbose = FALSE) {
  scanVCF(vcf, DP, GQ, samples, bannedPositions, regions = regions,
          verbose = verbose, returnGenotypeMatrix = FALSE)$callrate
}

#' Scan VCF file for genotypes.
#' 
#' The method scans VCF files and returns genotype matrix for specified 
#' variants after applying neccessary filters.
#' @inheritParams scanVCF
#' @param ... other arguments to be passed to scanVCF function.
#' @examples
#' vcf <- "CEU.exon.2010_09.genotypes.vcf.gz"
#' filepath <- system.file("extdata", vcf, package = "SVDFunctions")
#' 
#' samples <- c("NA06989", "NA10847", "NA11840", "NA12873")
#' genotypeMatrixVCF(filepath, GQ = 0, samples = samples) 
#' @export
genotypeMatrixVCF <- function(vcf, DP = 10L, GQ = 20L, variants = NULL,
                              samples = NULL, bannedPositions = NULL,
                              predictMissing = FALSE, 
                              verbose = FALSE, ...) {
  scanVCF(vcf, DP = DP, GQ = GQ, samples = samples, 
          bannedPositions = bannedPositions, variants = variants,
          verbose = verbose, predictMissing = predictMissing, ...)$genotype
}

#' Scan VCF file for sample names
#' 
#' The method scans only header of VCF file and returns a vector with sample 
#' names
#' @inheritParams scanVCF
#' @return character vector with sample names
#' @examples
#' vcf <- "CEU.exon.2010_09.genotypes.vcf.gz"
#' filepath <- system.file("extdata", vcf, package = "SVDFunctions")
#' 
#' sampleNamesVCF(filepath)
#' @export
sampleNamesVCF <- function(vcf, verbose = FALSE) {
    scanVCF(vcf, returnGenotypeMatrix = FALSE, verbose = verbose)$sample
}

#' Scan VCF file.
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
#' @param samples the set of samples to be scanned and returned.
#' @param bannedPositions the set of positions in format "chr#:#" that 
#' must be eliminated from consideration. 
#' @param variants the set of variants in format "chr#:# REF ALT"
#' (i.e. chr23:1532 T GT). In case of deletion ALT must be "*".
#' @param returnGenotypeMatrix logical: whether or not return genotype matrix.
#' @param predictMissing if TRUE missing values will be replaced with predicted
#' values.
#' @param missingRateThreshold the rows of VCF file with missing rate more
#' than threshold will be filtered out.
#' @param regions the set of regions in format "chr# startPos endPos". For each
#' region call rate will be calculated and corresponding matrix will be returned. 
#' @param binaryPathPrefix the path prefix for binary file prefix_bin and 
#' metadata file prefix_meta. If not NULL corresponding files will be generated.
#' @param verbose logical 
#' @return list containing genotype matrix and/or call rate matrix if 
#' requested.
#' @export
scanVCF <- function(vcf, DP = 10L, GQ = 20L, samples = NULL,
                    bannedPositions = NULL, variants = NULL, 
                    returnGenotypeMatrix = TRUE, predictMissing = FALSE, 
                    missingRateThreshold = 0.1, 
                    regions = NULL, binaryPathPrefix = NULL,
                    verbose = FALSE, ntree = 12) {
  vcf <- normalizePath(vcf)
  stopifnot(length(DP) > 0)
  stopifnot(length(GQ) > 0)
  DP <- as.integer(DP)
  GQ <- as.integer(GQ)
  stopifnot(!is.na(DP[0]))
  stopifnot(!is.na(GQ[0]))
  
  stopifnot(file.exists(vcf))
  
  tbi <- paste0(vcf, ".tbi")
  
  if (!is.null(binaryPathPrefix) | !file.exists(tbi)) {
    tbi <- NULL
  } else {
    if ((returnGenotypeMatrix && length(variants) != 0 && !predictMissing) |
        length(regions) != 0) {
      vcf <- createVCFFromTabixIndex(vcf, variants, regions, verbose)
    } else {
      tbi <- NULL
    }
  }
  
  fixChar <- function(x) if(is.null(x)) character(0) else x
  samples <- fixChar(samples)
  bannedPositions <- fixChar(bannedPositions)
  variants <- fixChar(variants)
  regions <- fixChar(regions)
  binaryPathPrefix <- fixChar(binaryPathPrefix)
  
  tryCatch( 
    res <- parse_vcf(vcf, samples, bannedPositions, variants, DP, GQ, 
                     returnGenotypeMatrix, isTRUE(predictMissing), regions, 
                     binaryPathPrefix, missingRateThreshold, ntree),
    error = function(c) {
      suffix <- ""
      if (!is.null(tbi)) {
        suffix <- paste0("Note: since .tbi is provided a temp file containing ",
                "only necessary regions has been created(", vcf, "). See ",
                "the temp file for more information about the error.\n")
      }
      stop(paste0(conditionMessage(c), "\n", suffix))
  })
  
  if (verbose) {
    print(data.frame(stat = names(res$stats), value = unlist(res$stats), 
                     row.names = NULL))
  }
  
  if (!is.null(res$genotype)) {
      colnames(res$genotype$genotype) <- res$samples
  }
  if (!isTRUE(predictMissing)) {
      res$genotype <- res$genotype$genotype
  }
  if (!is.null(res$callrate)) {
      colnames(res$callrate) <- res$samples
  }
  res
}

#' Collect counts statistics from genotype matrix.
#' 
#' Calculate counts statistics for genotype matrix and return data.frame
#' @param genotypeMatrix matrix containing values 0 for sample with both 
#' reference alleles, 1 for heterozygous samples, and 2 - with both alternative
#' alleles. Missing values are allowed.
#' @return data.frame with counts statistics.
#' @export
genotypesToCounts <- function(genotypeMatrix) {
  df <- NULL
  for (i in 0:2) {
    df <- cbind(df, apply(genotypeMatrix, 1, function(x) {
      sum(x[!is.na(x)] == i)
    }))
  }
  if (nrow(genotypeMatrix) == 0 || ncol(genotypeMatrix) == 0) {
    df <- matrix(ncol = 3, nrow = 0)
  }
  colnames(df) <- c("hom_ref", "het", "hom_alt")
  if (!is.null(rownames(genotypeMatrix))) {
    rownames(df) <- rownames(genotypeMatrix)
  }
  df
}

#' Scan binary variant file.
#' 
#' Scan binary file produced by \code{\link{scanVCF}} and collect required 
#' per variant allele counts statistics after applying quality control filters.
#' @inheritParams scanVCF
#' @param binaryFile the name of binary file.
#' @param metafile the name of metadata file.
#' @param samples the set of samples to be analyzed.
#' @param regions the set of regions [startPos, endPos] in format 
#' "chr# startPos endPos". Regions must be non-overlapping.
#' @param minMAF numeric minimum minor allele frequency.
#' @param maxMAF numeric maximum minor allele frequency.
#' @param minMAC integer minimum minor allele count.
#' @param maxMAC integer maximum minor allele count.
#' @param minCallRate numeric minimum call rate.
#' @param reportSingletons logical if TRUE singletons will be reported.
#' @return matrix with three columns:  
#' number of samples with 
#' both reference alleles, with  one reference allele and one alternative allele,
#' and with both alternative ones respectively.
#' @export
scanBinaryFile <- function(binaryFile, metafile, samples, 
                           variants = NULL, regions = NULL, 
                           DP = 10, GQ = 20, 
                           minMAF = 0.0, maxMAF = 1.0, 
                           minMAC = 1L, maxMAC = .Machine$integer.max, 
                           reportSingletons = TRUE, 
                           minCallRate = 0.9) {
  
  binaryFile <- normalizePath(binaryFile)
  metafile <- normalizePath(metafile)
  stopifnot(file.exists(binaryFile))
  stopifnot(file.exists(metafile))
  minMAC <- as.integer(minMAC)
  maxMAC <- as.integer(maxMAC)
  variants <- if (is.null(variants)) character(0) else variants  
  regions <- if (is.null(regions)) character(0) else regions
  
  res <- parse_binary_file(variants, samples, regions, binaryFile, metafile, 
                           minMAF, maxMAF, minCallRate, minMAC, maxMAC, 
                           reportSingletons,  DP, GQ);
  names <- res[["names"]]
  total <- as.integer(res["total"])
  res[["names"]] <- NULL
  res[["total"]] <- NULL
  res <- matrix(do.call(c, res), ncol = 4, dimnames = list(NULL, names(res)))
  alleles <- res[, "hom_ref"] + res[, "het"] + res[, "hom_alt"]
  res <- cbind(res, call_rate = alleles / (res[, "n_variants"] * total))
  rownames(res) <- names
  
  res
}
