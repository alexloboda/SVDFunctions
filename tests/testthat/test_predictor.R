context("replacing missing values") 

vcf <- "CEU.exon.2010_09.genotypes.vcf.gz"
file <- system.file("extdata", vcf, package = "SVDFunctions")

checkpoint_manifest_path <- function(dir) {
  file.path(dir, "predict_missing_checkpoint_predict_state.tsv")
}

checkpoint_previous_manifest_path <- function(dir) {
  file.path(dir, "predict_missing_checkpoint_predict_state.previous.tsv")
}

with_checkpoint_env <- function(vars, code) {
  old <- stats::setNames(as.list(Sys.getenv(names(vars), unset = NA_character_)), names(vars))
  on.exit({
    for (name in names(vars)) {
      if (is.na(old[[name]])) {
        Sys.unsetenv(name)
      } else {
        do.call(Sys.setenv, stats::setNames(list(old[[name]]), name))
      }
    }
  }, add = TRUE)

  for (name in names(vars)) {
    do.call(Sys.setenv, stats::setNames(list(vars[[name]]), name))
  }
  force(code)
}

checkpoint_predict_args <- function() {
  samples <- sampleNamesVCF(file)[1:20]
  base_gt <- genotypeMatrixVCF(file, DP = 20, GQ = 0, samples = samples)
  vars <- unique(c("chr1:76650391\tC\tT", rownames(base_gt)[1:min(12, nrow(base_gt))]))
  list(
    vcf = file,
    DP = 20,
    GQ = 0,
    samples = samples,
    variants = vars,
    predictMissing = TRUE,
    seed = 123L,
    window_size = 5L,
    rf_ntrees = 9L,
    checkpointInterval = 2L
  )
}

run_checkpoint_predict <- function(checkpointDir = NULL, resumeCheckpoint = TRUE, ...) {
  args <- utils::modifyList(checkpoint_predict_args(), list(...))
  if (!is.null(checkpointDir)) {
    args$checkpointDir <- checkpointDir
    args$resumeCheckpoint <- resumeCheckpoint
  }
  do.call(genotypeMatrixVCF, args)
}
  
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

test_that("predictMissing checkpoint resumes after simulated crash during segment write", {
  expected <- run_checkpoint_predict()
  checkpoint_dir <- tempfile("predict-checkpoint-")
  dir.create(checkpoint_dir)

  expect_error(
    with_checkpoint_env(
      c(SVDF_CHECKPOINT_TEST_ABORT_AFTER_SEGMENT_WRITE = "2"),
      run_checkpoint_predict(checkpointDir = checkpoint_dir)
    ),
    "Simulated checkpoint crash after segment write"
  )

  resumed <- run_checkpoint_predict(checkpointDir = checkpoint_dir)
  expect_equal(resumed$genotype, expected$genotype)
  expect_equal(resumed$predicted, expected$predicted)
  expect_equal(resumed$loo, expected$loo)
})

test_that("predictMissing checkpoint falls back to previous manifest when current is corrupt", {
  checkpoint_dir <- tempfile("predict-checkpoint-")
  dir.create(checkpoint_dir)

  expected <- run_checkpoint_predict(checkpointDir = checkpoint_dir)
  current_manifest <- checkpoint_manifest_path(checkpoint_dir)
  previous_manifest <- checkpoint_previous_manifest_path(checkpoint_dir)
  expect_true(file.exists(current_manifest))
  expect_true(file.exists(previous_manifest))

  writeLines("broken-current-manifest", current_manifest)

  resumed <- expect_warning(
    run_checkpoint_predict(checkpointDir = checkpoint_dir),
    "recovered using"
  )
  expect_equal(resumed$genotype, expected$genotype)
  expect_equal(resumed$predicted, expected$predicted)
  expect_equal(resumed$loo, expected$loo)
})

test_that("predictMissing checkpoint errors when the only visible checkpoint is invalid", {
  checkpoint_dir <- tempfile("predict-checkpoint-")
  dir.create(checkpoint_dir)

  expect_error(
    with_checkpoint_env(
      c(SVDF_CHECKPOINT_TEST_ABORT_AFTER_SAVE = "1"),
      run_checkpoint_predict(checkpointDir = checkpoint_dir)
    ),
    "Simulated checkpoint crash after manifest save"
  )

  current_manifest <- checkpoint_manifest_path(checkpoint_dir)
  manifest_lines <- readLines(current_manifest)
  manifest_lines <- sub("^resume_offset\t.*$", "resume_offset\t999999999999", manifest_lines)
  writeLines(manifest_lines, current_manifest)

  expect_error(
    run_checkpoint_predict(checkpointDir = checkpoint_dir),
    "Checkpoint resume offset|seek BGZF stream|past the last data row"
  )
})

test_that("predictMissing complete checkpoint is reused without rewriting state", {
  checkpoint_dir <- tempfile("predict-checkpoint-")
  dir.create(checkpoint_dir)

  expected <- run_checkpoint_predict(checkpointDir = checkpoint_dir)
  current_manifest <- checkpoint_manifest_path(checkpoint_dir)
  manifest_before <- paste(readLines(current_manifest), collapse = "\n")

  reused <- with_checkpoint_env(
    c(SVDF_CHECKPOINT_TEST_ABORT_AFTER_SEGMENT_WRITE = "1"),
    run_checkpoint_predict(checkpointDir = checkpoint_dir)
  )
  manifest_after <- paste(readLines(current_manifest), collapse = "\n")

  expect_equal(reused$genotype, expected$genotype)
  expect_equal(reused$predicted, expected$predicted)
  expect_equal(reused$loo, expected$loo)
  expect_identical(manifest_after, manifest_before)
})
