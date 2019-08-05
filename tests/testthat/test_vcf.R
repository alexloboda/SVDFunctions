context("parsing VCF files")

vcf <- "CEU.exon.2010_09.genotypes.vcf.gz"
file <- system.file("extdata", vcf, package = "SVDFunctions")
DP <- 20

missing_function <- function(x) {
  startsWith(x, ".") | 
  sapply(strsplit(x, ":"), function(x) as.integer(x[2]) < DP)
}

test_that("genotype matrix are parsed correctly", {
  bannedPos <- "chr1:18022097"
  banned <- "chr1:18022097\tG\tT"
  samples <- sampleNamesVCF(file)
  vcf <- genotypeMatrixVCF(file, DP = DP, GQ = 0, samples = samples[1:10], 
                               bannedPositions = bannedPos)
  df <- seqminer::tabix.read.table(file, paste0(c(1:22), ":1-300000000"))
  format <- function(x) paste0("chr", df$CHROM[x], ":", df$POS[x], "\t", 
                               df$REF[x], "\t", df$ALT[x])
  rownames(df) <- lapply(1:nrow(df), format)
  df <- df[, samples[1:10]]
  df <- apply(df, c(1, 2), function(x) if(missing_function(x)) NA else x)
  df <- apply(df, c(1, 2), function(x) {
    if (is.na(x)) {
      NA
    } else if (startsWith(x, "0/1") | startsWith(x, "1/0")) {
      1
    } else if (startsWith(x, "0/0")) {
      0
    } else if (startsWith(x, "1/1")) {
      2
    }
  })
  df <- df[rowSums(is.na(df)) <= 1, ]
  df <- df[rowSums(matrix(!is.na(df) & df == 1, nrow = nrow(df))) > 0 | 
          (rowSums(matrix(!is.na(df) & df == 0, nrow = nrow(df))) > 0 & 
           rowSums(matrix(!is.na(df) & df == 2, nrow = nrow(df))) > 0), ]
  df <- df[rownames(df) != banned, ]
  expect_equal(df, vcf)
})

test_that("callrates are calculated correctly", {
  regions <- data.frame(chr = c("1", "1", "20"), 
                        from = c("1108138", "40000000", "33521213"), 
                        to = c("36000000", "40000010", "61665700"))
  samples <- c("NA06989", "NA10847", "NA11840", "NA12873")
  pkgFormat = function(x) paste0("chr", x[1], " ", x[2], " ", x[3])
  seqMinerFormat = function(x) paste0(x[1], ":", x[2], "-", x[3])
  regionsPkg <- apply(regions, 1, pkgFormat)
  regionsSeqMiner <- apply(regions, 1, seqMinerFormat)
  cr <- callRateMatrixVCF(file, DP = DP, GQ = 0, regions = regionsPkg, 
                           samples = samples)
  expectedMatrix <- matrix(nrow = 0, ncol = length(samples))
  for (i in 1:nrow(regions)) {
    df <- seqminer::tabix.read.table(file, regionsSeqMiner[i])
    if (nrow(df) == 0) {
      expectedMatrix <- rbind(expectedMatrix, NA)
      next
    }
    df <- df[, samples]
    missing <- sapply(df[, samples], function(x) sum(missing_function(x)))
    expected <- (nrow(df) - missing) / nrow(df)
    expectedMatrix <- rbind(expectedMatrix, expected)
  }
  rownames(expectedMatrix) <- regionsPkg
  expectedMatrix <- expectedMatrix[apply(expectedMatrix, 1, 
                                         function(x) !any(is.na(x))), ]
  expect_equal(cr, expectedMatrix, tolerance = 1e-8)
})

test_that("storing/extracting data to/from binary file works ", {
  set.seed(42)
  prefix <- paste0(tempdir(), "/db")
  bin <- paste0(prefix, "_bin")
  meta <- paste0(prefix, "_meta")
  vcf <- scanVCF(file, DP = DP, GQ = 0, binaryPath = prefix)
  
  samples <- vcf$samples
  variants <- rownames(vcf$genotype)
  variants <- c(variants, "chr1:1\nT\nC")
  samples <- sample(samples, as.integer(length(samples) / 2))
  sample_size <- length(variants) %/% 2
  sel <- c(rep(TRUE, sample_size), rep(FALSE, length(variants) - sample_size))
  variants <- variants[sample(sel)]
  
  localDP <- 30
  vcf <- scanVCF(file, DP = localDP, GQ = 0, samples = samples, variants = variants)
  actual <- scanBinaryFile(bin, meta, samples, variants, DP = localDP, GQ = 0)
  expected <- genotypesToCounts(vcf$genotype)
  expect_equal(actual[, 1:3], expected)
  
  regions <- c("chr2 148943492 234267016",
               "chr7 102949810 149848267",
               "chr10 46506985 46507378")
  
  actualReg <- list()
  actualReg[["complete"]] <- t(matrix(c(c(96, 42, 27), 
                                        c(84, 35, 4), 
                                        c(38, 45, 0)), nrow = 3))
  
  actualReg[["maf"]] <- t(matrix(c(c(17, 39, 27), 
                                   c(44, 34, 4), 
                                   c(38, 45, 0)), nrow = 3))
  
  actualReg[["rare"]] <- t(matrix(c(c(79, 3, 0),
                                    c(40, 1, 0),
                                    c(0, 0, 0)), nrow = 3))
  
  actualReg[["cr"]] <- t(matrix(c(c(17, 22, 3), 
                                  c(0, 0, 0),
                                  c(18, 24, 0)), nrow = 3))
  
  actualReg <- lapply(actualReg, function(x) {
    rownames(x) <- regions
    colnames(x) <- c("hom_ref", "het", "hom_alt")
    x <- x[apply(x, 1, function(x) sum(x) > 0), ]
    x
  })
  
  reg <- list()
  reg[["complete"]] <- scanBinaryFile(bin, meta, samples, regions = regions, 
                        DP = localDP, GQ = 0)[, 1:3]
  reg[["maf"]] <- scanBinaryFile(bin, meta, samples, regions = regions, 
                        DP = localDP, GQ = 0, minMAF = 0.04)[, 1:3]
  reg[["rare"]] <- scanBinaryFile(bin, meta, samples, regions = regions, 
                        DP = localDP, GQ = 0, maxMAF = 0.04)[, 1:3]
  reg[["cr"]] <- scanBinaryFile(bin, meta, samples, regions = regions, 
                        DP = localDP, GQ = 0, minCallRate = 41.5 / 45)[, 1:3]
  
  expect_equal(actualReg, reg)
})

test_that("parsing multivariant lines works", {
  file <- system.file("extdata", "multivariant.vcf.gz", package = "SVDFunctions")
  GT <- genotypeMatrixVCF(file, DP = 0, GQ = 0)
  expected <- matrix(c(0, 1, 1, 1, 1, 2), ncol = 3)
  colnames(expected) <- c("A", "B", "C")
  rownames(expected) <- c("chr1:1\tT\tG", "chr1:2\tT\t*")
  expect_equal(GT, expected)
})
