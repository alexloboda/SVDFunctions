#' Intermediate computation of the residual vectors 
#' 
#' Residual vector is estimated as \eqn{(I - UU^T)Z}, where part 
#' \eqn{I - UU^T} is the same for every 
#' control sample. Thus, could be precomputed only once for every case 
#' basis supplied. 
#' @param SV Number of singular vectors to be used for reconstruction
#' of the original vector
#' @param referenceU Matrix of the left singular vectors of cases
#' @export
computeResidual.preproc <- function (referenceU, SV) 
{
    return(diag(length(referenceU[, 1])) - tcrossprod(x = referenceU[, 
        SV]))
}
