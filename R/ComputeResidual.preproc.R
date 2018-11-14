#' Intermediate computation of the residual vectors 
#' 
#' Residual vector is estimated as \eqn{(I - UU^T)Z}, where part 
#' \eqn{I - UU^T} is the same for every 
#' control sample. Thus, could be pre-computed only once for every case 
#' basis supplied. 
#' @param SV Number of singular vectors to be used for reconstruction
#' of the original vector
#' @param ReferenceU Matrix of the left singular vectors of cases
#' @export
ComputeResidual.preproc <- function (ReferenceU, SV) 
{
    return(diag(length(ReferenceU[, 1])) - tcrossprod(x = ReferenceU[, 
        SV]))
}
