# Moore-Penrose pseudo-inverse computed from the singular value decomposition.
# This is a dependency-free, numerically stable replacement for pracma::pinv
# and uses the same default tolerance for discarding singular values.
pinv <- function(m, tol = NULL) {
  s <- svd(m)
  if (is.null(tol)) {
    tol <- max(dim(m)) * .Machine$double.eps * max(s$d)
  }
  positive <- s$d > tol
  if (!any(positive)) {
    return(matrix(0, nrow = ncol(m), ncol = nrow(m)))
  }
  s$v[, positive, drop = FALSE] %*%
    ((1 / s$d[positive]) * t(s$u[, positive, drop = FALSE]))
}

# Seeds for the C++ side. Drawing the default from R's own generator is what
# makes set.seed() reach the simulated annealing: the C++ code has no access to
# R's RNG state, so the seed has to be handed over explicitly.
resolveSeed <- function(seed) {
  if (is.null(seed)) {
    seed <- stats::runif(1, 1, .Machine$integer.max)
  }
  seed <- as.integer(seed)
  stopifnot(length(seed) == 1, !is.na(seed), seed > 0)
  seed
}

# Column names used inside ggplot2::aes() in the plotting helpers. Declaring
# them keeps R CMD check from reporting "no visible binding" notes.
utils::globalVariables(c(
  "PC1", "PC2", "startx", "starty", "endx", "endy",
  "expected", "cupper", "clower", "X"
))

