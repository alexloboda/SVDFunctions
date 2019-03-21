#' Estimation of residual vector norms for all controls
#' 
#' @param genotypeMatrix Genotype matrix
#' @param SVDReference Reference basis of the left singular vectors
#' @param nSV Number of singular vectors to be used for reconstruction of the 
#' original vector
#' @export
parallelResidEstimate <- function (genotypeMatrix, SVDReference, nSV) 
{
    gmatrix <- genotypeMatrix
    preprocessed.ref <- ComputeResidual.preproc(SVDReference, seq(1, nSV, 1))
    resiudals <- preprocessed.ref %*% gmatrix
    apply(resiudals, MARGIN = 2, function(x) norm(x, type = "2"))
}
