context("file backed control selection")

# The same construction the other selectControls tests use, plus a raw genotype
# matrix that actually carries missing values: the file backed path has to keep
# the imputed and the raw matrix apart exactly as the in-memory one does.
buildFileFixture <- function(nclust = 0L) {
  set.seed(20240607)

  p <- 300L
  n_controls <- 160L
  n_cases <- 60L
  components <- 8L

  af <- runif(p, 0.15, 0.45)
  drawGenotypes <- function(n) {
    t(vapply(seq_len(p), function(i) rbinom(n, 2, af[i]), integer(n)))
  }

  controls <- drawGenotypes(n_controls)
  storage.mode(controls) <- "integer"
  rownames(controls) <- paste0("v", seq_len(p))
  colnames(controls) <- paste0("ctrl", seq_len(n_controls))

  cases <- drawGenotypes(n_cases)
  rownames(cases) <- rownames(controls)
  colnames(cases) <- paste0("case", seq_len(n_cases))

  # genotypeMatrix stays complete, as selectControls requires; the raw matrix
  # loses a few genotypes, which is what the counting stage has to cope with.
  genotypeMatrix <- controls
  storage.mode(genotypeMatrix) <- "double"
  raw <- controls
  missing <- sample(length(raw), 200L)
  raw[missing] <- NA_integer_

  controlsMean <- rowMeans(controls)
  svdRef <- RSpectra::svds(controls - controlsMean, k = components)$u
  rownames(svdRef) <- rownames(controls)

  reducedCases <- t(svdRef) %*% (cases - controlsMean)
  casesMean <- rowMeans(reducedCases)
  caseSvd <- svd(reducedCases - casesMean)
  casesPDs <- (caseSvd$u %*% diag(caseSvd$d)) / sqrt(n_cases)

  clustering <- if (nclust > 0L) {
    paste0("cl", rep_len(seq_len(nclust), n_controls))
  } else {
    NULL
  }

  list(genotypeMatrix = genotypeMatrix, raw = raw, casesPDs = casesPDs,
       casesMean = casesMean, svdRef = svdRef, controlsMean = controlsMean,
       caseCounts = genotypesToCounts(cases), clustering = clustering,
       variants = rownames(controls), samples = colnames(controls))
}

writeFixture <- function(fx, dir, dtype = "f64") {
  paths <- list(genotype = file.path(dir, "genotype.ctlm"),
                raw = file.path(dir, "raw.ctlm"),
                collapsed = file.path(dir, "collapsed.ctlm"))
  writeControlsMatrix(fx$genotypeMatrix, paths$genotype, dtype = dtype)
  writeControlsMatrix(fx$raw, paths$raw, dtype = "u8")
  collapseControlsMatrix(paths$raw, paths$collapsed,
                         controlsClustering = fx$clustering)
  paths
}

test_that("a control matrix survives a write/read round trip", {
  fx <- buildFileFixture()
  dir <- tempfile("ctlm"); dir.create(dir)
  path <- file.path(dir, "genotype.ctlm")
  writeControlsMatrix(fx$genotypeMatrix, path, dtype = "f64")

  info <- readControlsMatrixInfo(path)
  expect_equal(info$dtype, "f64")
  expect_equal(info$nrow, nrow(fx$genotypeMatrix))
  expect_equal(info$ncol, ncol(fx$genotypeMatrix))
  expect_identical(info$rownames, fx$variants)
  expect_identical(info$colnames, fx$samples)
  # Every column of a genotype matrix stands for exactly one sample.
  expect_identical(info$colweights, rep(1L, ncol(fx$genotypeMatrix)))

  rawPath <- file.path(dir, "raw.ctlm")
  writeControlsMatrix(fx$raw, rawPath, dtype = "u8")
  expect_equal(readControlsMatrixInfo(rawPath)$dtype, "u8")

  fractional <- fx$genotypeMatrix
  fractional[1, 1] <- 0.5
  expect_error(writeControlsMatrix(fractional, file.path(dir, "bad.ctlm"),
                                   dtype = "u8"),
               "genotypes 0, 1, 2 or NA")
  noNames <- fx$genotypeMatrix
  dimnames(noNames) <- NULL
  expect_error(writeControlsMatrix(noNames, file.path(dir, "bad.ctlm")),
               "variant row names")
})

test_that("the streamed projection matches the in-memory one", {
  fx <- buildFileFixture()
  dir <- tempfile("ctlm"); dir.create(dir)
  path <- file.path(dir, "genotype.ctlm")
  writeControlsMatrix(fx$genotypeMatrix, path, dtype = "f64")

  inMemory <- prepareControlsSpace(fx$genotypeMatrix, fx$svdRef, fx$controlsMean,
                                   fx$casesMean, fx$casesPDs)
  streamed <- prepareControlsSpace(path, fx$svdRef, fx$controlsMean,
                                   fx$casesMean, fx$casesPDs)

  # The two paths sum the same products in a different order, so they agree to
  # rounding rather than bit for bit.
  expect_equal(streamed$points, inMemory$points, tolerance = 1e-12)
  expect_lt(max(abs(streamed$points - inMemory$points)) / max(abs(inMemory$points)),
            1e-13)
  for (field in setdiff(names(inMemory), "points")) {
    expect_identical(streamed[[field]], inMemory[[field]], info = field)
  }
})

test_that("the streamed projection does not depend on the thread count", {
  fx <- buildFileFixture()
  dir <- tempfile("ctlm"); dir.create(dir)
  path <- file.path(dir, "genotype.ctlm")
  writeControlsMatrix(fx$genotypeMatrix, path, dtype = "f64")

  project <- function(threads) {
    prepareControlsSpace(path, fx$svdRef, fx$controlsMean, fx$casesMean,
                         fx$casesPDs, threads = threads)$points
  }
  reference <- project(1L)
  for (threads in c(2L, 3L, 4L, 8L)) {
    expect_identical(project(threads), reference)
  }
})

test_that("the streamed projection handles a variant subset in any row order", {
  fx <- buildFileFixture()
  dir <- tempfile("ctlm"); dir.create(dir)

  shuffled <- fx$genotypeMatrix[sample(nrow(fx$genotypeMatrix)), ]
  path <- file.path(dir, "shuffled.ctlm")
  writeControlsMatrix(shuffled, path, dtype = "f64")

  caseVariants <- fx$variants[seq_len(200L)]
  inMemory <- prepareControlsSpace(shuffled, fx$svdRef, fx$controlsMean,
                                   fx$casesMean, fx$casesPDs,
                                   caseVariants = caseVariants)
  streamed <- prepareControlsSpace(path, fx$svdRef, fx$controlsMean,
                                   fx$casesMean, fx$casesPDs,
                                   caseVariants = caseVariants)

  expect_equal(streamed$points, inMemory$points, tolerance = 1e-12)
  expect_identical(streamed$variantRows, inMemory$variantRows)
})

test_that("the projection refuses a matrix with missing genotypes", {
  fx <- buildFileFixture()
  dir <- tempfile("ctlm"); dir.create(dir)
  path <- file.path(dir, "raw.ctlm")
  writeControlsMatrix(fx$raw, path, dtype = "u8")

  expect_error(prepareControlsSpace(path, fx$svdRef, fx$controlsMean,
                                    fx$casesMean, fx$casesPDs),
               "must not be missing")

  collapsed <- file.path(dir, "collapsed.ctlm")
  collapseControlsMatrix(path, collapsed)
  expect_error(prepareControlsSpace(collapsed, fx$svdRef, fx$controlsMean,
                                    fx$casesMean, fx$casesPDs),
               "not collapsed cluster counts")
})

test_that("matching is identical from a matrix, a raw file and a collapsed file", {
  for (nclust in c(0L, 40L)) {
    fx <- buildFileFixture(nclust)
    dir <- tempfile("ctlm"); dir.create(dir)
    paths <- writeFixture(fx, dir)

    space <- prepareControlsSpace(fx$genotypeMatrix, fx$svdRef, fx$controlsMean,
                                  fx$casesMean, fx$casesPDs,
                                  controlsClustering = fx$clustering)
    sizes <- if (nclust > 0L) list(12L, 20L, 2L) else list(60L, 100L, 10L)
    candidates <- mvnSubsampleClusters(space$points, space$mean, space$cov,
                                       space$clusterIds, min = sizes[[1]],
                                       max = sizes[[2]], step = sizes[[3]],
                                       iterations = 2000L, seed = 4242L,
                                       threads = 1L, precomputeThreads = 1L)

    args <- list(caseCounts = fx$caseCounts, candidates = candidates,
                 minLambda = 0.5, softMinLambda = 0.9, softMaxLambda = 1.05,
                 maxLambda = 1.3, min = sizes[[1]], minCallRate = 0.9)
    fromMatrix <- do.call(matchControlCandidates,
                          c(list(fx$raw, clusterIds = space$clusterIds), args))
    fromRaw <- do.call(matchControlCandidates,
                       c(list(paths$raw, clusterIds = space$clusterIds), args))
    fromCollapsed <- do.call(matchControlCandidates,
                             c(list(paths$collapsed, clusterIds = space$clusterIds),
                               args))

    # Counting is exact arithmetic, so the three sources must agree bit for bit.
    expect_identical(fromRaw, fromMatrix, info = paste("nclust", nclust))
    expect_identical(fromCollapsed, fromMatrix, info = paste("nclust", nclust))
    expect_true(length(fromMatrix$controls) > 0)
  }
})

test_that("a collapsed file carries its own clustering", {
  fx <- buildFileFixture(40L)
  dir <- tempfile("ctlm"); dir.create(dir)
  paths <- writeFixture(fx, dir)

  info <- readControlsMatrixInfo(paths$collapsed)
  expect_equal(info$dtype, "clusterCounts")
  expect_equal(info$ncol, 40L)
  expect_identical(info$colnames, sort(unique(fx$clustering)))
  expect_equal(sum(info$colweights), ncol(fx$raw))

  space <- prepareControlsSpace(fx$genotypeMatrix, fx$svdRef, fx$controlsMean,
                                fx$casesMean, fx$casesPDs,
                                controlsClustering = fx$clustering)
  # Same labels, in the same order, as the projection stage resolved.
  expect_identical(info$colnames, space$clusterLabels)

  candidates <- mvnSubsampleClusters(space$points, space$mean, space$cov,
                                     space$clusterIds, min = 12L, max = 20L,
                                     step = 2L, iterations = 2000L, seed = 7L,
                                     threads = 1L, precomputeThreads = 1L)
  supplied <- matchControlCandidates(paths$collapsed, fx$caseCounts, candidates,
                                     clusterIds = space$clusterIds,
                                     minLambda = 0.5, min = 12L, minCallRate = 0.9)
  derived <- matchControlCandidates(paths$collapsed, fx$caseCounts, candidates,
                                    minLambda = 0.5, min = 12L, minCallRate = 0.9)
  expect_identical(derived, supplied)
})

# The streamed and the in-memory projection agree only to rounding, and the
# annealing that consumes them is a discrete search: a last-bit difference in the
# coordinates sends it down a different, equally valid trajectory. So the wiring
# is pinned down exactly against a pipeline that shares the same coordinates, and
# the comparison against selectControls is a comparison of quality.
test_that("selectControlsFromFiles wires the three stages together", {
  resolveSelectedControls <- getFromNamespace("resolveSelectedControls",
                                              "SVDFunctions")

  for (nclust in c(0L, 40L)) {
    fx <- buildFileFixture(nclust)
    dir <- tempfile("ctlm"); dir.create(dir)
    paths <- writeFixture(fx, dir)
    sizes <- if (nclust > 0L) list(12L, 20L, 2L) else list(60L, 100L, 10L)
    label <- paste("nclust", nclust)

    common <- list(casesPDs = fx$casesPDs, casesMean = fx$casesMean,
                   SVDReference = fx$svdRef, controlsMean = fx$controlsMean,
                   caseCounts = fx$caseCounts, controlsClustering = fx$clustering,
                   minLambda = 0.5, softMinLambda = 0.9, softMaxLambda = 1.05,
                   maxLambda = 1.3, min = sizes[[1]], max = sizes[[2]],
                   step = sizes[[3]], iterations = 2000L, minCallRate = 0.9,
                   seed = 4242L, saThreads = 1L, exactPrecomputeThreads = 1L)

    fromRaw <- do.call(selectControlsFromFiles,
                       c(list(genotypeFile = paths$genotype,
                              originalFile = paths$raw), common))

    # Same stage one (same file, so the same coordinates), stages two and three
    # run by hand against the in-memory matrix.
    space <- prepareControlsSpace(paths$genotype, fx$svdRef, fx$controlsMean,
                                  fx$casesMean, fx$casesPDs,
                                  controlsClustering = fx$clustering)
    candidates <- mvnSubsampleClusters(space$points, space$mean, space$cov,
                                       space$clusterIds, min = sizes[[1]],
                                       max = sizes[[2]], step = sizes[[3]],
                                       iterations = 2000L, seed = 4242L,
                                       threads = 1L, precomputeThreads = 1L)
    byHand <- matchControlCandidates(fx$raw, fx$caseCounts, candidates,
                                     space$clusterIds, space$variantRows,
                                     minLambda = 0.5, softMinLambda = 0.9,
                                     softMaxLambda = 1.05, maxLambda = 1.3,
                                     min = sizes[[1]], minCallRate = 0.9)
    byHand$controls <- resolveSelectedControls(byHand$controls,
                                               space$clusterIds - 1L,
                                               space$clusterLabels,
                                               colnames(fx$raw), FALSE)
    expect_identical(fromRaw, byHand, info = label)
    expect_true(length(fromRaw$controls) > 0)

    # A collapsed file selects the very same clusters, and reports them as such.
    fromCollapsed <- do.call(selectControlsFromFiles,
                             c(list(genotypeFile = paths$genotype,
                                    originalFile = paths$collapsed), common))
    asClusters <- do.call(selectControlsFromFiles,
                          c(list(genotypeFile = paths$genotype,
                                 originalFile = paths$raw,
                                 returnClusters = TRUE), common))
    expect_identical(fromCollapsed, asClusters, info = label)
  }
})

test_that("the file backed selection is as good as the in-memory one", {
  fx <- buildFileFixture(40L)
  dir <- tempfile("ctlm"); dir.create(dir)
  paths <- writeFixture(fx, dir)

  common <- list(casesPDs = fx$casesPDs, casesMean = fx$casesMean,
                 SVDReference = fx$svdRef, controlsMean = fx$controlsMean,
                 caseCounts = fx$caseCounts, controlsClustering = fx$clustering,
                 minLambda = 0.5, softMinLambda = 0.9, softMaxLambda = 1.05,
                 maxLambda = 1.3, min = 12L, max = 20L, step = 2L,
                 iterations = 5000L, minCallRate = 0.9, seed = 4242L,
                 saThreads = 1L, exactPrecomputeThreads = 1L)

  inMemory <- do.call(selectControls,
                      c(list(genotypeMatrix = fx$genotypeMatrix,
                             originalGenotypeMatrix = fx$raw), common))
  fromFiles <- do.call(selectControlsFromFiles,
                       c(list(genotypeFile = paths$genotype,
                              originalFile = paths$raw), common))

  expect_true(length(inMemory$controls) > 0)
  expect_equal(length(fromFiles$controls), length(inMemory$controls))
  expect_true(all(fromFiles$controls %in% colnames(fx$raw)))
  expect_lt(as.numeric(fromFiles$optimal_lambda)[1], 1.3)
  expect_gt(as.numeric(fromFiles$optimal_lambda)[1], 0.5)
  # Both searches see the same coordinates up to rounding, so they should land on
  # comparably good solutions even when the subsets differ.
  expect_equal(as.numeric(fromFiles$optimal_lambda)[1],
               as.numeric(inMemory$optimal_lambda)[1], tolerance = 0.1)
})

test_that("f32 storage stays within single precision of the exact projection", {
  fx <- buildFileFixture()
  dir <- tempfile("ctlm"); dir.create(dir)
  path <- file.path(dir, "genotype32.ctlm")
  writeControlsMatrix(fx$genotypeMatrix, path, dtype = "f32")
  expect_equal(readControlsMatrixInfo(path)$dtype, "f32")

  inMemory <- prepareControlsSpace(fx$genotypeMatrix, fx$svdRef, fx$controlsMean,
                                   fx$casesMean, fx$casesPDs)
  streamed <- prepareControlsSpace(path, fx$svdRef, fx$controlsMean,
                                   fx$casesMean, fx$casesPDs)
  relative <- max(abs(streamed$points - inMemory$points)) / max(abs(inMemory$points))
  expect_lt(relative, 1e-6)
})

test_that("mismatched files are rejected instead of silently misaligned", {
  fx <- buildFileFixture(40L)
  dir <- tempfile("ctlm"); dir.create(dir)
  paths <- writeFixture(fx, dir)

  common <- list(casesPDs = fx$casesPDs, casesMean = fx$casesMean,
                 SVDReference = fx$svdRef, controlsMean = fx$controlsMean,
                 caseCounts = fx$caseCounts, minLambda = 0.5, min = 12L,
                 max = 20L, step = 2L, iterations = 200L, minCallRate = 0.9,
                 seed = 1L)

  # A different number of clusters changes the labels the file stores.
  expect_error(do.call(selectControlsFromFiles,
                       c(list(genotypeFile = paths$genotype,
                              originalFile = paths$collapsed,
                              controlsClustering =
                                paste0("cl", rep_len(seq_len(20), ncol(fx$raw)))),
                         common)),
               "different")

  # Same labels, different sizes: caught by the sample counts the file carries.
  lopsided <- c(rep("cl1", 10L),
                rep_len(paste0("cl", 2:40), ncol(fx$raw) - 10L))
  expect_error(do.call(selectControlsFromFiles,
                       c(list(genotypeFile = paths$genotype,
                              originalFile = paths$collapsed,
                              controlsClustering = lopsided), common)),
               "different")

  # Reordered samples in the raw file must not pass as the same controls.
  swapped <- fx$raw[, c(2:ncol(fx$raw), 1L)]
  swappedPath <- file.path(dir, "swapped.ctlm")
  writeControlsMatrix(swapped, swappedPath, dtype = "u8")
  expect_error(do.call(selectControlsFromFiles,
                       c(list(genotypeFile = paths$genotype,
                              originalFile = swappedPath,
                              controlsClustering = fx$clustering), common)),
               "same samples in the same order")

  expect_error(collapseControlsMatrix(paths$genotype, file.path(dir, "x.ctlm")),
               "dtype \"u8\"")
})

test_that("collapsing rejects a cluster that does not fit in a byte", {
  dir <- tempfile("ctlm"); dir.create(dir)
  wide <- matrix(0L, nrow = 2L, ncol = 300L,
                 dimnames = list(c("v1", "v2"), paste0("s", seq_len(300L))))
  path <- file.path(dir, "wide.ctlm")
  writeControlsMatrix(wide, path, dtype = "u8")

  expect_error(collapseControlsMatrix(path, file.path(dir, "out.ctlm"),
                                      controlsClustering = rep("one", 300L)),
               "at most 255 samples")
  expect_silent(collapseControlsMatrix(path, file.path(dir, "out.ctlm"),
                                       controlsClustering =
                                         rep(c("a", "b"), each = 150L)))
})
