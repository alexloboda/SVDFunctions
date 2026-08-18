context("hierarchical control selection")

# Two leaves and the split that merges them: the smallest instance that runs
# every branch of recSelect and mergedOrJoint.
buildHierFixture <- function(nclust = 0L, uneven = FALSE) {
  set.seed(31337)
  p <- 300L; n_controls <- 160L; n_case <- 45L; components <- 8L

  af <- runif(p, 0.15, 0.45)
  draw <- function(n) t(vapply(seq_len(p), function(i) rbinom(n, 2, af[i]), integer(n)))

  controls <- draw(n_controls)
  storage.mode(controls) <- "integer"
  rownames(controls) <- paste0("v", seq_len(p))
  colnames(controls) <- paste0("ctrl", seq_len(n_controls))

  controlsMean <- rowMeans(controls)
  svdRef <- RSpectra::svds(controls - controlsMean, k = components)$u
  rownames(svdRef) <- rownames(controls)
  transition <- t(svdRef)

  makePopulation <- function(g, id) {
    reduced <- transition %*% (g - controlsMean)
    m <- rowMeans(reduced)
    s <- svd(reduced - m)
    list(US = (s$u %*% diag(s$d)) / sqrt(ncol(g)), mean = m,
         counts = genotypesToCounts(g), id = id, cluster = id)
  }

  a <- draw(n_case); rownames(a) <- rownames(controls)
  b <- draw(n_case); rownames(b) <- rownames(controls)

  cases <- list(
    variants = rownames(controls),
    population = list(makePopulation(a, 1L), makePopulation(b, 2L),
                      makePopulation(cbind(a, b), 3L)),
    hierarchy = list(left = list(id = 1L, name = "A", type = "leaf"),
                     right = list(id = 2L, name = "B", type = "leaf"),
                     type = "split", id = 3L),
    names = stats::setNames(c("A", "B", "A & B"), c("1", "2", "3")))

  gm <- controls; storage.mode(gm) <- "double"
  clustering <- if (uneven) {
    # Deliberately unequal clusters, so weighting them by size is observable.
    rep(paste0("cl", seq_len(20L)), times = rep(c(4L, 12L), 10L))
  } else if (nclust > 0L) {
    paste0("cl", rep_len(seq_len(nclust), n_controls))
  } else {
    NULL
  }

  list(genotypeMatrix = gm, controls = controls, cases = cases, svdRef = svdRef,
       controlsMean = controlsMean, clustering = clustering)
}

hierArgs <- function(fx, nclust) {
  sizes <- if (nclust > 0L) list(6L, 12L, 2L) else list(60L, 100L, 10L)
  list(cases = fx$cases, SVDReference = fx$svdRef, controlsMean = fx$controlsMean,
       controlsClustering = fx$clustering, minLambda = 0.5, softMinLambda = 0.9,
       softMaxLambda = 1.05, maxLambda = 1.3, min = sizes[[1]], max = sizes[[2]],
       step = sizes[[3]], iterations = 2000L, minCallRate = 0.9, seed = 4242L,
       saThreads = 1L, exactPrecomputeThreads = 1L)
}

test_that("retargetControlsSpace keeps the coordinates and swaps the target", {
  fx <- buildHierFixture()
  space <- prepareControlsSpace(fx$genotypeMatrix, fx$svdRef, fx$controlsMean,
                                fx$cases$population[[1]]$mean,
                                fx$cases$population[[1]]$US)
  second <- fx$cases$population[[2]]
  retargeted <- retargetControlsSpace(space, second$mean, second$US)

  expect_identical(retargeted$points, space$points)
  expect_identical(retargeted$clusterIds, space$clusterIds)
  expect_identical(retargeted$variantRows, space$variantRows)
  expect_identical(retargeted$mean, second$mean)
  expect_identical(retargeted$cov, tcrossprod(second$US))

  # Building the space around that population directly must give the same thing.
  direct <- prepareControlsSpace(fx$genotypeMatrix, fx$svdRef, fx$controlsMean,
                                 second$mean, second$US)
  expect_identical(retargeted, direct)

  expect_error(retargetControlsSpace(space, second$mean[1:3], second$US),
               "must match the components")
  expect_error(retargetControlsSpace(list(), second$mean, second$US),
               "as returned by prepareControlsSpace")
})

test_that("prepared cluster counts match the one-shot path and can be reused", {
  fx <- buildHierFixture(40L)
  space <- prepareControlsSpace(fx$genotypeMatrix, fx$svdRef, fx$controlsMean,
                                fx$cases$population[[1]]$mean,
                                fx$cases$population[[1]]$US,
                                controlsClustering = fx$clustering)
  candidates <- mvnSubsampleClusters(space$points, space$mean, space$cov,
                                     space$clusterIds, min = 12L, max = 20L,
                                     step = 2L, iterations = 2000L, seed = 11L,
                                     threads = 1L, precomputeThreads = 1L)

  counts <- prepareClusterCounts(fx$controls, space$clusterIds, space$variantRows)
  expect_s3_class(counts, "clusterCounts")
  expect_equal(counts$nVariants, length(space$variantRows))
  expect_equal(counts$nClusters, 40L)
  expect_output(print(counts), "clusterCounts")

  args <- list(candidates = candidates, minLambda = 0.5, min = 12L,
               minCallRate = 0.9)
  # The same counts scored against each population must equal the one-shot call.
  for (population in fx$cases$population) {
    prepared <- do.call(matchControlCandidates,
                        c(list(counts, population$counts), args))
    oneShot <- do.call(matchControlCandidates,
                       c(list(fx$controls, population$counts,
                              clusterIds = space$clusterIds,
                              variantRows = space$variantRows), args))
    expect_identical(prepared, oneShot)
  }

  expect_error(matchControlCandidates(counts, fx$cases$population[[1]]$counts[1:5, ],
                                      candidates),
               "one row per case variant")
})

test_that("the hierarchy prepares the control side once for every node", {
  fx <- buildHierFixture()
  ns <- asNamespace("SVDFunctions")
  calls <- new.env()
  calls$space <- 0L
  calls$counts <- 0L
  suppressMessages({
    trace("prepareControlsSpace", where = ns, print = FALSE,
          tracer = function() calls$space <- calls$space + 1L)
    trace("prepareClusterCounts", where = ns, print = FALSE,
          tracer = function() calls$counts <- calls$counts + 1L)
  })
  on.exit(suppressMessages({
    untrace("prepareControlsSpace", where = ns)
    untrace("prepareClusterCounts", where = ns)
  }), add = TRUE)

  do.call(selectControlsHier,
          c(list(fx$genotypeMatrix, fx$controls), hierArgs(fx, 0L)))

  expect_equal(length(fx$cases$population), 3L)
  expect_equal(calls$space, 1L)
  expect_equal(calls$counts, 1L)
})

test_that("the merge heuristic does not depend on how controls are reported", {
  fx <- buildHierFixture(uneven = TRUE)
  # Wide soft bounds so at least one node counts as good and ncontrols is
  # actually exercised rather than short-circuited to zero.
  args <- utils::modifyList(hierArgs(fx, 20L),
                            list(softMinLambda = 0.5, softMaxLambda = 1.5))

  bySample <- do.call(selectControlsHier,
                      c(list(fx$genotypeMatrix, fx$controls), args))
  byCluster <- do.call(selectControlsHier,
                       c(list(fx$genotypeMatrix, fx$controls),
                         args, list(returnClusters = TRUE)))

  # Clusters here hold 4 or 12 samples, so counting cluster identifiers instead
  # of weighing them by size would give a different -- and wrong -- answer.
  expect_identical(byCluster$ncontrols, bySample$ncontrols)
  expect_identical(byCluster$clusters, bySample$clusters)
  expect_identical(byCluster$lambdas, bySample$lambdas)
  expect_identical(byCluster$minL, bySample$minL)
  expect_gt(bySample$ncontrols, 0)
  expect_lt(length(unique(byCluster$table$sample)),
            length(unique(bySample$table$sample)))
})

test_that("selectControlsHierFromFiles matches the in-memory hierarchy", {
  for (nclust in c(0L, 20L)) {
    fx <- buildHierFixture(nclust)
    dir <- tempfile("ctlm"); dir.create(dir)
    genotypeFile <- file.path(dir, "genotype.ctlm")
    rawFile <- file.path(dir, "raw.ctlm")
    collapsedFile <- file.path(dir, "collapsed.ctlm")
    writeControlsMatrix(fx$genotypeMatrix, genotypeFile, dtype = "f64")
    writeControlsMatrix(fx$controls, rawFile, dtype = "u8")
    collapseControlsMatrix(rawFile, collapsedFile,
                           controlsClustering = fx$clustering)

    args <- hierArgs(fx, nclust)
    fromRaw <- do.call(selectControlsHierFromFiles,
                       c(list(genotypeFile, rawFile), args))
    inMemory <- do.call(selectControlsHier,
                        c(list(fx$genotypeMatrix, fx$controls), args))

    # Same coordinates up to rounding, so compare what has to agree exactly
    # against a run that shares them, and the rest by quality.
    asClusters <- do.call(selectControlsHierFromFiles,
                          c(list(genotypeFile, rawFile), args,
                            list(returnClusters = TRUE)))
    fromCollapsed <- do.call(selectControlsHierFromFiles,
                             c(list(genotypeFile, collapsedFile), args))
    expect_identical(fromCollapsed$table, asClusters$table)
    expect_identical(fromCollapsed$lambdas, asClusters$lambdas)
    expect_identical(fromCollapsed$ncontrols, asClusters$ncontrols)

    expect_equal(nrow(fromRaw$table), nrow(inMemory$table))
    expect_identical(names(fromRaw$lambdas), names(inMemory$lambdas))
    expect_identical(fromRaw$names, inMemory$names)
    expect_true(all(fromRaw$table$sample %in% colnames(fx$controls)))
    expect_lt(fromRaw$minL, 1.3)
  }
})

test_that("hierarchy file inputs are validated", {
  fx <- buildHierFixture(20L)
  dir <- tempfile("ctlm"); dir.create(dir)
  genotypeFile <- file.path(dir, "genotype.ctlm")
  rawFile <- file.path(dir, "raw.ctlm")
  collapsedFile <- file.path(dir, "collapsed.ctlm")
  writeControlsMatrix(fx$genotypeMatrix, genotypeFile, dtype = "f64")
  writeControlsMatrix(fx$controls, rawFile, dtype = "u8")
  collapseControlsMatrix(rawFile, collapsedFile,
                         controlsClustering = fx$clustering)

  args <- hierArgs(fx, 20L)
  args$controlsClustering <- paste0("x", rep_len(seq_len(10L),
                                                 ncol(fx$controls)))
  expect_error(do.call(selectControlsHierFromFiles,
                       c(list(genotypeFile, collapsedFile), args)),
               "different")

  args$controlsClustering <- fx$clustering
  args$cases$variants <- paste0("zz", seq_along(args$cases$variants))
  expect_error(do.call(selectControlsHierFromFiles,
                       c(list(genotypeFile, rawFile), args)),
               "subset of the control")
})
