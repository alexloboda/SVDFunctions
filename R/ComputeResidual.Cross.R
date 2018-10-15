#' Computing relative residual vector norm
#' 
#' Computes relative residual vector norm using pre-computed results 
#' from ComputeResidual.preproc
#' @param SampleVector Vector of genotypes for a sample to be evaluated as control
#' @param  PreprocessedU \eqn{I - UU^T} matrix calculated with ComputeResidual.preproc
#' @export
ComputeResidual.Cross <- function (SampleVector, PreprocessedU) 
{
    Residual <- tcrossprod(PreprocessedU, t(SampleVector))
    Residual <- norm(Residual, type = "2")/norm(SampleVector, 
        type = "2")
    return(Residual)
}

