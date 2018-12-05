context("constructing an instance from raw data")

test_that("necessary files are created without errors", {
  rawDataPath <- system.file("extdata", package = "SVDFunctions")
  tmp <- tempdir()
  bfilename <- "regions_extracted"
  bfile <- paste0(rawDataPath, "/", bfilename)
  ref <- paste0(rawDataPath, "/list_of_variants_available_for_matching_with_ref_allele.txt")
  BED2GMatrix(bfile, ref, outputDir = tmp)
  
  U <- data.table::fread(paste0(tmp, "/", bfilename, "_U.txt"))
  gmatrix <- data.table::fread(paste0(tmp, "/", bfilename, "_gmatrix.txt"))
  case_counts <- data.table::fread(paste0(tmp, "/", bfilename, "_case_counts.txt"))
  
  expect_equal(dim(U), c(4028, 10))
  expect_equal(dim(gmatrix), c(4028, 49))
  expect_equal(dim(case_counts), c(4028, 6))
})