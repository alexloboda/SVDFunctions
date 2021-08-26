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
#' @export
normal_subsample <- function(matrix, size) {
  subsample_mvn(matrix, size, colMeans(matrix), stats::cov(t(matrix)))
}