context("parsing VCF files")

vcf <- "CEU.exon.2010_09.genotypes.vcf.gz"

test_that("callrates are calculated correctly", {
  file <- system.file("extdata", vcf, package = "SVDFunctions")
  DP <- 20
  regions <- data.frame(chr = c("1", "20"), 
                        from = c("1108138", "33521213"), 
                        to = c("36000000", "61665700"))
  pkgFormat = function(x) paste0("chr", x[1], " ", x[2], " ", x[3])
  seqMinerFormat = function(x) paste0(x[1], ":", x[2], "-", x[3])
  regionsPkg <- apply(regions, 1, pkgFormat)
  regionsSeqMiner <- apply(regions, 1, seqMinerFormat)
  samples <- c("NA06989", "NA10847", "NA11840", "NA12873")
  vcf <- scanVCF(file, DP = DP, GQ = 0, regions = regionsPkg, samples = samples)
  expectedMatrix <- matrix(nrow = 0, ncol = length(samples))
  for (i in 1:nrow(regions)) {
    df <- seqminer::tabix.read.table(file, regionsSeqMiner[i])[, samples]
    missing_function <- function(x) {
      startsWith(x, ".") | 
      sapply(strsplit(x, ":"), function(x) as.integer(x[2]) < DP)
    }
    missing <- sapply(df[, samples], function(x) sum(missing_function(x)))
    expected <- (nrow(df) - missing) / nrow(df)
    expectedMatrix <- rbind(expectedMatrix, expected)
  }
  rownames(expectedMatrix) <- regionsPkg
  expect_equal(vcf$callrate, expectedMatrix, tolerance = 1e-8)
})
