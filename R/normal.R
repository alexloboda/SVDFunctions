#' Find a subset of points forming multivariate normal distribution.
#' 
#' The method takes a set of n-dimensional points of size m and tries 
#' to select \code{size} points from the set by minimizing difference
#' between emperical characteristic function of subset and characteristic funciton
#' of multivariate normal with sample mean and covariance from the whole set
#' as parameters.
#' @param matrix n by m R matrix where columns are samples and rows are
#' corresponding components of vectors.
#' @param size the number of points to be subsetted.
#' @param seed optional integer seed for the simulated annealing. When
#' \code{NULL} (default) the seed is drawn from R's generator, so
#' \code{set.seed} makes the result reproducible. Reproducibility holds for a
#' given build of the package on a given machine; see \code{\link{selectControls}}.
#' @export
normal_subsample <- function(matrix, size, seed = NULL) {
  subsample_mvn(matrix, size, rowMeans(matrix), stats::cov(t(matrix)),
                resolveSeed(seed))
}