context("replacing missing values") 

vcf <- "CEU.exon.2010_09.genotypes.vcf.gz"
file <- system.file("extdata", vcf, package = "SVDFunctions")
  
test_that("predictMissing flag works on one variant from 1kG", {
  gm <- genotypeMatrixVCF(file, DP = 20, GQ = 0, predictMissing = TRUE,
                          variants = "chr1:76650391\tC\tT")
  expect_false(any(is.na(gm$genotype)))
  expect_equal(sum(gm$predicted), 6)
})

test_that("predictMissing flag works on all variants from 1kG", {
  gm <- genotypeMatrixVCF(file, DP = 20, GQ = 0)
  gm <- genotypeMatrixVCF(file, DP = 20, GQ = 0, predictMissing = TRUE)
  expect_false(any(is.na(gm$genotype)))
})

test_that("predictMissing is deterministic across repeated calls", {
  samples <- sampleNamesVCF(file)[1:20]
  base_gt <- genotypeMatrixVCF(file, DP = 20, GQ = 0, samples = samples)
  vars <- rownames(base_gt)[1:min(5, nrow(base_gt))]

  gm_one <- genotypeMatrixVCF(
    file,
    DP = 20,
    GQ = 0,
    samples = samples,
    variants = vars,
    predictMissing = TRUE
  )
  gm_two <- genotypeMatrixVCF(
    file,
    DP = 20,
    GQ = 0,
    samples = samples,
    variants = vars,
    predictMissing = TRUE
  )

  expect_equal(gm_one$genotype, gm_two$genotype)
  expect_equal(gm_one$predicted, gm_two$predicted)
  expect_equal(gm_one$loo, gm_two$loo)
})

test_that("predictMissing accepts rf_ntrees and keeps default behavior", {
  samples <- sampleNamesVCF(file)[1:20]
  base_gt <- genotypeMatrixVCF(file, DP = 20, GQ = 0, samples = samples)
  vars <- rownames(base_gt)[1:min(5, nrow(base_gt))]

  gm_default <- genotypeMatrixVCF(
    file,
    DP = 20,
    GQ = 0,
    samples = samples,
    variants = vars,
    predictMissing = TRUE
  )
  gm_explicit <- genotypeMatrixVCF(
    file,
    DP = 20,
    GQ = 0,
    samples = samples,
    variants = vars,
    predictMissing = TRUE,
    rf_ntrees = 50L
  )
  gm_small <- genotypeMatrixVCF(
    file,
    DP = 20,
    GQ = 0,
    samples = samples,
    variants = vars,
    predictMissing = TRUE,
    rf_ntrees = 5L
  )

  expect_equal(gm_default$genotype, gm_explicit$genotype)
  expect_equal(gm_default$predicted, gm_explicit$predicted)
  expect_equal(gm_default$loo, gm_explicit$loo)

  expect_equal(dim(gm_small$genotype), dim(gm_default$genotype))
  expect_equal(dim(gm_small$predicted), dim(gm_default$predicted))
  expect_equal(colnames(gm_small$loo), colnames(gm_default$loo))
})
