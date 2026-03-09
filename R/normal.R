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
  stopifnot(is.matrix(matrix))
  size <- as.integer(size)

  if (nrow(matrix) == 0 || ncol(matrix) == 0) {
    stop("Matrix must be non-empty.")
  }
  if (length(size) != 1 || is.na(size)) {
    stop("size must be a single integer value.")
  }
  if (size < nrow(matrix) + 1 || size > ncol(matrix)) {
    stop("Requested subsample size must be between nrow(matrix) + 1 and ncol(matrix).")
  }

  subsample_mvn(matrix, size, rowMeans(matrix), stats::cov(t(matrix)))
}