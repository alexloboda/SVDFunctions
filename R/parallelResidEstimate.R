#' Estimation of residual vector norms for all controls.
#' 
#' @param genotypeMatrix Genotype matrix.
#' @param SVDReference Reference basis of the left singular vectors.
#' @export
parallelResidEstimate <- function(genotypeMatrix, SVDReference,
                                   caseMeans,  contMeans, SV = 3) {
  SVint <- ncol(SVDReference)
  caseMeans <- caseMeans - contMeans
  genotypeMatrix <- genotypeMatrix - contMeans
  contDecomp <- RSpectra::svds(genotypeMatrix, k = SV)
  basisChange <- t(contDecomp$u)
  center_in_un <- basisChange %*% matrix(caseMeans, ncol = 1)
  covMatrix <- basisChange %*% SVDReference %*%
    t(SVDReference) %*% t(basisChange)
  casesCovInv <- solve(covMatrix, tol = 1e-30)
  mahalanobis <- function(vector) {
    new_vector <- matrix(basisChange %*% vector - center_in_un, ncol = 1)
    t(new_vector) %*% (casesCovInv %*% new_vector)
  }
  scores <- apply(genotypeMatrix, 2, mahalanobis)
  scores
}
