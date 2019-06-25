#' Estimation of residual vector norms for all controls
#' 
#' @param genotypeMatrix Genotype matrix
#' @param SVDReference Reference basis of the left singular vectors
#' @param nSV Number of singular vectors to be used for reconstruction of the 
#' original vector
#' @export
parallelResidEstimate <- function (genotypeMatrix, SVDReference, nSV) 
{
  SVDReference <- SVDReference[, 1:nSV]
  inter <- crossprod(SVDReference, genotypeMatrix) 
  ret <- c()
  idx <- 1:ncol(genotypeMatrix)
  groups <- split(idx, ceiling(seq_along(idx) / 256))
  for (group in groups) {
    start <- group[1] - 1
    rhs <- SVDReference %*% inter[, group]
    for (i in group) {
      v <- genotypeMatrix[, i] - rhs[, i - start]
      ret <- c(ret, sqrt(sum(v * v)))
      rm(v)
    }
    rm(rhs)
  }
  names(ret) <- colnames(genotypeMatrix)
  ret
}
