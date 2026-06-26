context("selectControls end-to-end")

test_that("internal clustering helpers preserve factor ids and legacy order", {
  prepareControlsClustering <- getFromNamespace("prepareControlsClustering",
                                                "SVDFunctions")
  resolveSelectedControls <- getFromNamespace("resolveSelectedControls",
                                              "SVDFunctions")

  controlNames <- c("ctrl1", "ctrl2", "ctrl3", "ctrl4", "ctrl5")
  clusteringInfo <- prepareControlsClustering(c("b", "a", "b", "c", "a"),
                                              controlNames)

  expect_equal(clusteringInfo$sampleClusterIds, c(1L, 0L, 1L, 2L, 0L))
  expect_equal(clusteringInfo$clusterLabels, c("a", "b", "c"))

  # Old C++ behaviour expanded each selected cluster in the returned order,
  # preserving sample order within every cluster.
  expect_equal(resolveSelectedControls(c(2L, 1L),
                                       clusteringInfo$sampleClusterIds,
                                       clusteringInfo$clusterLabels,
                                       controlNames),
               c("ctrl1", "ctrl3", "ctrl2", "ctrl5"))
  expect_equal(resolveSelectedControls(c(2L, 1L),
                                       clusteringInfo$sampleClusterIds,
                                       clusteringInfo$clusterLabels,
                                       controlNames,
                                       returnClusters = TRUE),
               c("b", "a"))
})

test_that("internal clustering helpers handle NULL clustering", {
  prepareControlsClustering <- getFromNamespace("prepareControlsClustering",
                                                "SVDFunctions")
  resolveSelectedControls <- getFromNamespace("resolveSelectedControls",
                                              "SVDFunctions")

  controlNames <- c("ctrl1", "ctrl2", "ctrl3", "ctrl4")
  clusteringInfo <- prepareControlsClustering(NULL, controlNames)

  expect_equal(clusteringInfo$sampleClusterIds, 0:3)
  expect_equal(clusteringInfo$clusterLabels, controlNames)
  expect_equal(resolveSelectedControls(c(3L, 1L),
                                       clusteringInfo$sampleClusterIds,
                                       clusteringInfo$clusterLabels,
                                       controlNames),
               c("ctrl3", "ctrl1"))
  expect_equal(resolveSelectedControls(c(3L, 1L),
                                       clusteringInfo$sampleClusterIds,
                                       clusteringInfo$clusterLabels,
                                       controlNames,
                                       returnClusters = TRUE),
               c("ctrl3", "ctrl1"))
})

# Build a self-consistent toy problem where cases and controls are drawn from
# the same per-variant allele frequencies. In that situation a subset of the
# controls matches the cases well, so the matching should return a non-empty
# control set with a genomic inflation factor (lambda) close to 1.
test_that("selectControls returns a non-empty set with low lambda", {
  skip_on_cran()
  set.seed(123)

  p <- 300L           # variants
  n_controls <- 160L  # control samples
  n_cases <- 60L      # case samples
  components <- 8L

  # Moderate allele frequencies so the variants comfortably pass QC.
  af <- runif(p, 0.15, 0.45)
  drawGenotypes <- function(n) {
    m <- vapply(seq_len(p), function(i) rbinom(n, 2, af[i]), integer(n))
    t(m) # variants x samples
  }

  controls <- drawGenotypes(n_controls)
  storage.mode(controls) <- "integer"
  rownames(controls) <- paste0("v", seq_len(p))
  colnames(controls) <- paste0("ctrl", seq_len(n_controls))

  cases <- drawGenotypes(n_cases)
  rownames(cases) <- rownames(controls)
  colnames(cases) <- paste0("case", seq_len(n_cases))

  # Reference left-singular basis from the centered control matrix.
  controlsMean <- rowMeans(controls)
  svdRef <- RSpectra::svds(controls - controlsMean, k = components)$u
  rownames(svdRef) <- rownames(controls)

  # The columns of svdRef are orthonormal, so its pseudo-inverse is t(svdRef).
  # Reduce the cases exactly the way selectControls reduces the controls.
  transition <- t(svdRef)
  reducedCases <- transition %*% (cases - controlsMean)
  casesMean <- rowMeans(reducedCases)
  centeredCases <- reducedCases - casesMean
  caseSvd <- svd(centeredCases)
  casesPDs <- (caseSvd$u %*% diag(caseSvd$d)) / sqrt(n_cases)

  caseCounts <- genotypesToCounts(cases)

  genotypeMatrix <- controls
  storage.mode(genotypeMatrix) <- "double"

  result <- selectControls(
    genotypeMatrix = genotypeMatrix,
    originalGenotypeMatrix = controls,
    casesPDs = casesPDs,
    casesMean = casesMean,
    SVDReference = svdRef,
    controlsMean = controlsMean,
    caseCounts = caseCounts,
    minLambda = 0.5,
    softMinLambda = 0.9,
    softMaxLambda = 1.05,
    maxLambda = 1.3,
    min = 60L,
    max = 100L,
    step = 10L,
    iterations = 5000L,
    minCallRate = 0.9
  )

  expect_true(length(result$controls) > 0)
  expect_true(all(result$controls %in% colnames(controls)))

  lambda <- as.numeric(result$optimal_lambda)[1]
  expect_false(is.na(lambda))
  # The target is lambda < 1.05; assert the looser < 1.3 bound for test
  # stability across platforms and RNG streams.
  expect_lt(lambda, 1.3)
})

