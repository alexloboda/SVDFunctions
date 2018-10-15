#' Estimation of relative residual vector norms for all controls
#' 
#' @param gmatrix Genotype matrix
#' @param svdReference Reference basis of the left singular vectors
#' @param nSV Number of singular vectors to be used for reconstruction of the 
#' original vector
#' @export
ParallelResidEstimate <- function (gmatrix, svdReference, nSV) 
{
    preprocessed.ref <- ComputeResidual.preproc(svdReference, seq(1, nSV, 1))
    resiudals <- preprocessed.ref %*% gmatrix
    norms <- c()
    for (i in 1:ncol(resiudals)) {
        normRes <- norm(resiudals[, i], type = "2")
        normG <- norm(gmatrix[, i], type = "2")
        norms <- append(norms, normRes)
    }
    norms
}
