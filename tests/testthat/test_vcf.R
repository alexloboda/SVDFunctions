context("parsing VCF files")


test_that("parsing gives reasonable genotype matrix", {
  file <- system.file("extdata", "CEU.exon.2010_09.genotypes.vcf.gz",
                      package = "SVDFunctions")
  samples <- c("NA07051", "NA12045", "NA12400")
  variants <- c("chr1:3537996\tT\tC", "chr1:3538692\tG\tC", "chr5:141290008\tG\tA",
                "chr11:64281842\tC\tT")
  banned <- c("chr1:3538692")
  expected <- t(matrix(c(2, 2, NA, 2, 1, 0, 1, 0, NA), nrow = length(samples), 
                     dimnames = list(samples, variants[c(1,3,4)])))
  vcf <- scanVCF(file, DP = 20, GQ = 0, samples = samples, 
                 bannedPositions = banned, variants = variants)
  expect_equal(vcf$genotype, expected)
})

test_that("callrates are calculated correctly", {
  file <- system.file("extdata", "CEU.exon.2010_09.genotypes.vcf.gz",
                      package = "SVDFunctions")
  regions <- c("chr1 1108138 36000000")
  samples <- c("NA07051", "NA12045", "NA12400")
  vcf <- scanVCF(file, DP = 20, GQ = 0, regions = regions, 
               samples = samples)
  expectedCallrate <- matrix(c(1.0, 1.0, 3.0 / 14.0), nrow = length(regions), 
                            dimnames = list(regions, samples))
  expect_equal(vcf$callrate, expectedCallrate, tolerance = 1e-8)
})
