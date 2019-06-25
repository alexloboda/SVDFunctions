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
  residuals <- genotypeMatrix - rhs
  ret <- c()
  for (i in 1:ncol(residuals)) {
    ret <- c(ret, norm(residuals[, i], type = "2"))
  }
  names(ret) <- colnames(residuals)
  ret
}
