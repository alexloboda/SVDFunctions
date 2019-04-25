context("replacing missing values") 

vcf <- "CEU.exon.2010_09.genotypes.vcf.gz"
file <- system.file("extdata", vcf, package = "SVDFunctions")
  
test_that("predictMissing flag works on one variant from 1kG", {
  gm <- genotypeMatrixVCF(file, DP = 20, GQ = 0, predictMissing = TRUE,
                          variants = "chr1:76650391\tC\tT")
  expect_false(any(is.na(gm)))
})

test_that("predictMissing flag works on all variants from 1kG", {
  gm <- genotypeMatrixVCF(file, DP = 20, GQ = 0)
  gm <- genotypeMatrixVCF(file, DP = 20, GQ = 0, predictMissing = TRUE, 
                          variants = rownames(gm)[1:50])
  expect_false(any(is.na(gm)))
})

test_that("excluding samples from predictor works", {
  vcf <- "predictor_test_dataset.vcf.gz"
  file <- system.file("extdata", vcf, package = "SVDFunctions")
  gm <- genotypeMatrixVCF(file, DP = 0, GQ = 0, predictMissing = TRUE)
  expect_equal(gm[1, 2], 1.0, tolerance = 0.2)
  
  gm <- genotypeMatrixVCF(file, DP = 0, GQ = 0, predictMissing = TRUE, 
                          excludedPredictorSamples = c("F", "G", "H", "I", "J"))
  expect_equal(gm[1, 2], 1.57, tolerance = 0.1)
})
