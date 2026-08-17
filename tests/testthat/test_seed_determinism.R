context("selectControls seed determinism")

# Builds a small but non-degenerate control/case setup, the same construction the
# lambda test uses, and returns everything selectControls needs. Using real PCA
# projections (rather than a toy Gaussian) keeps the normality statistic positive
# so the annealing performs real moves instead of a random walk.
buildFixture <- function() {
  set.seed(20240607)

  p <- 300L           # variants
  n_controls <- 160L  # control samples
  n_cases <- 60L      # case samples
  components <- 8L

  af <- runif(p, 0.15, 0.45)
  drawGenotypes <- function(n) {
    m <- vapply(seq_len(p), function(i) rbinom(n, 2, af[i]), integer(n))
    t(m)
  }

  controls <- drawGenotypes(n_controls)
  storage.mode(controls) <- "integer"
  rownames(controls) <- paste0("v", seq_len(p))
  colnames(controls) <- paste0("ctrl", seq_len(n_controls))

  cases <- drawGenotypes(n_cases)
  rownames(cases) <- rownames(controls)
  colnames(cases) <- paste0("case", seq_len(n_cases))

  controlsMean <- rowMeans(controls)
  svdRef <- RSpectra::svds(controls - controlsMean, k = components)$u
  rownames(svdRef) <- rownames(controls)

  transition <- t(svdRef)
  reducedCases <- transition %*% (cases - controlsMean)
  casesMean <- rowMeans(reducedCases)
  centeredCases <- reducedCases - casesMean
  caseSvd <- svd(centeredCases)
  casesPDs <- (caseSvd$u %*% diag(caseSvd$d)) / sqrt(n_cases)

  genotypeMatrix <- controls
  storage.mode(genotypeMatrix) <- "double"

  list(genotypeMatrix = genotypeMatrix, controls = controls,
       casesPDs = casesPDs, casesMean = casesMean, svdRef = svdRef,
       controlsMean = controlsMean, caseCounts = genotypesToCounts(cases))
}

# Runs selectControls against the fixture with an explicit seed and thread layout
# and returns the selected control names.
selectWith <- function(fx, seed, saThreads = 0L, exactPrecomputeThreads = 0L) {
  selectControls(
    genotypeMatrix = fx$genotypeMatrix,
    originalGenotypeMatrix = fx$controls,
    casesPDs = fx$casesPDs,
    casesMean = fx$casesMean,
    SVDReference = fx$svdRef,
    controlsMean = fx$controlsMean,
    caseCounts = fx$caseCounts,
    minLambda = 0.5,
    softMinLambda = 0.9,
    softMaxLambda = 1.05,
    maxLambda = 1.3,
    min = 60L,
    max = 100L,
    step = 10L,
    iterations = 5000L,
    minCallRate = 0.9,
    seed = seed,
    saThreads = saThreads,
    exactPrecomputeThreads = exactPrecomputeThreads
  )$controls
}

test_that("a fixed seed reproduces the same controls", {
  skip_on_cran()
  fx <- buildFixture()

  reference <- selectWith(fx, seed = 4242L)
  expect_true(length(reference) > 0)

  # Loop, because the historical failure mode was a ~1-in-a-few flake that a
  # single comparison would routinely miss.
  for (i in seq_len(5)) {
    expect_identical(selectWith(fx, seed = 4242L), reference)
  }
})

test_that("the result is independent of the thread layout", {
  skip_on_cran()
  fx <- buildFixture()

  reference <- selectWith(fx, seed = 4242L, saThreads = 1L,
                          exactPrecomputeThreads = 1L)

  expect_identical(selectWith(fx, seed = 4242L, saThreads = 4L,
                              exactPrecomputeThreads = 1L), reference)
  expect_identical(selectWith(fx, seed = 4242L, saThreads = 1L,
                              exactPrecomputeThreads = 4L), reference)
  expect_identical(selectWith(fx, seed = 4242L, saThreads = 4L,
                              exactPrecomputeThreads = 4L), reference)
})

test_that("different seeds explore different subsets", {
  skip_on_cran()
  fx <- buildFixture()

  reference <- selectWith(fx, seed = 4242L)
  others <- vapply(c(1L, 99L, 100000L, 2000000000L),
                   function(s) identical(selectWith(fx, seed = s), reference),
                   logical(1))
  # The seed must actually steer the search: at least one of several distinct
  # seeds has to yield a different selection.
  expect_false(all(others))
})
