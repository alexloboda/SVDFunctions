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
  rhs <- SVDReference %*% inter
  ret <- c()
  for (i in 1:ncol(genotypeMatrix)) {
    v <- genotypeMatrix[, i] - rhs[, i]
    ret <- c(ret, norm(v, type = "2"))
    rm(v)
  }
  names(ret) <- colnames(genotypeMatrix)
  ret
}
