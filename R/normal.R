#' @export
normal_subsample <- function(matrix, size) {
  subsample_mvn(matrix, size, colMeans(matrix), stats::cov(t(matrix)))
}
